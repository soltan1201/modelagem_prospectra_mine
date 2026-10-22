"""
Treino da CNN de favorabilidade de fosfato com os patches de exportar_patches.py
(plano em plano_cnn.md). PyTorch + timm; roda igual no servidor, na GPU local
e no Colab (o notebook só chama main() ou as funções).

Entrada   PATCHES_DIR/modelo_<A|B1>/patches_<LxCy>.npz  (X (N, C, 64, 64) float16)
          A  = 15 índices + swir_estimado   B1 = 9 log-bandas + 5 emissividades + swir_estimado
Divisão   campos 'conjunto' e 'fold' dos .npz (divisao.json): treino nos folds ≠ k,
          validação no fold k; 'teste' só com --teste (uma vez, no fim).

Modelos (--backbone)
  pequena               CNN pequena do zero (conv + SE), para 32 ou 64 px
  regnety_016           RegNetY (atenção de canal SE embutida)        ┐ timm,
  tf_efficientnetv2_s   EfficientNetV2-S                              │ pesos ImageNet,
  resnet50              ResNet-50                                     │ 1ª camada adaptada
  coatnet_0_rw_224      CoAtNet-0 (híbrida conv + atenção)            ┘ a C canais (in_chans)
  Backbones pré-treinados: primeiras CONGELAR épocas só com a cabeça treinando,
  depois tudo com lr menor.

Rótulos   y contínuo (rótulo suave: média de p_RF/p_GTB nos pseudo-rótulos, 1 nas
          ocorrências) e peso por amostra → BCE com alvo suave ponderada. Peso
          extra equilibra positivos (y ≥ 0,5) e negativos.
Normalização  por canal, estatísticas só do treino (corte p1–p99, média, desvio);
          swir_estimado fica 0/1; NaN → 0 depois de padronizar.
Aumento   rotações de 90°, espelhamento, ganho/deslocamento aleatório por canal
          (simula outra cena — contra a assinatura de cena).

Métricas (validação)
  auc_amostra   todos os patches, positivo = y ≥ 0,5 (pseudo-rótulos)
  auc_ponto     ocorrências reais do fold (mediana dos jitters por ponto) × negativos
  ap, brier
  Referência: GTB, AUC por ponto com célula de fora = 0,84 (v1).

Saídas (--saida/<modelo>_<backbone>_<patch>px_<data>/)
  config.json, metricas_folds.csv, oof_validacao.csv (predição fora do fold),
  fold<k>.pt (melhor época), teste.csv / metricas_teste.json (com --teste)

Uso
  python treinar_cnn.py --modelo A --backbone pequena --patch 32
  python treinar_cnn.py --modelo B1 --backbone regnety_016 --folds 0,1,2,3,4
  python treinar_cnn.py --modelo A --backbone regnety_016 --teste     # final
"""

import argparse
import json
import os
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import average_precision_score, roc_auc_score
from torch.utils.data import DataLoader, Dataset

PATCHES_DIR = Path(os.environ.get('PATCHES_DIR', '/srv/almacen/db_images/patches_cnn')).expanduser()
RUNS_DIR    = Path(os.environ.get('RUNS_DIR', '/srv/almacen/db_images/runs_cnn')).expanduser()
CANAL_BIN   = 'swir_estimado'     # canal binário: não padroniza nem recebe ruído
GTB_REF     = 0.84


# ═══════════════════════════════════════════════════════════════════════════
# DADOS
# ═══════════════════════════════════════════════════════════════════════════
META = ['id', 'y', 'peso', 'origem_amostra', 'conjunto', 'fold', 'celula', 'lon', 'lat',
        'p_rf', 'p_gtb', 'papel', 'codigo_obj', 'jitter']


def carregar(modelo, pasta=PATCHES_DIR):
    """Junta os .npz das células: X (N, C, H, W) float16, metadados e nomes dos canais."""
    arquivos = sorted((pasta / f'modelo_{modelo}').glob('patches_*.npz'))
    if not arquivos:
        raise SystemExit(f'Nenhum .npz em {pasta / f"modelo_{modelo}"} — rode exportar_patches.py')
    Xs, metas, canais = [], [], None
    for a in arquivos:
        z = np.load(a, allow_pickle=True)    # .npz gerados por exportar_patches (texto pode vir como objeto)
        Xs.append(z['X'])
        metas.append(pd.DataFrame({k: z[k] for k in META if k in z.files}))
        canais = [str(c) for c in z['canais']]
    meta = pd.concat(metas, ignore_index=True)
    meta['fold'] = meta['fold'].astype(int)
    meta['y'] = meta['y'].astype(np.float32)
    meta['peso'] = meta['peso'].astype(np.float32)
    return np.concatenate(Xs), meta, canais


def estatisticas(X, canais, idx_treino, max_patches=3000, seed=0):
    """Corte p1–p99, média e desvio por canal, só com patches de treino."""
    rng = np.random.default_rng(seed)
    sel = rng.choice(idx_treino, min(max_patches, len(idx_treino)), replace=False)
    amostra = X[sel].astype(np.float32)
    est = {}
    for c, nome in enumerate(canais):
        v = amostra[:, c].ravel()
        v = v[np.isfinite(v)]
        if nome == CANAL_BIN or not len(v):
            est[nome] = (0.0, 1.0, 0.0, 1.0)
            continue
        lo, hi = np.percentile(v, [1, 99])
        v = np.clip(v, lo, hi)
        est[nome] = (float(lo), float(hi), float(v.mean()), float(v.std() + 1e-6))
    return est


class Patches(Dataset):
    def __init__(self, X, meta, idx, canais, est, patch, aumentar):
        self.X, self.meta, self.idx = X, meta, np.asarray(idx)
        self.canais, self.patch, self.aumentar = canais, patch, aumentar
        a = np.array([est[c] for c in canais], np.float32)
        self.lo, self.hi, self.mu, self.sd = (a[:, k][:, None, None] for k in range(4))
        self.bin = np.array([c == CANAL_BIN for c in canais])

    def __len__(self):
        return len(self.idx)

    def __getitem__(self, i):
        j = self.idx[i]
        bruto = self.X[j].astype(np.float32)
        o = (bruto.shape[-1] - self.patch) // 2            # recorte central (ex.: 32 de 64 px)
        bruto = bruto[:, o:o + self.patch, o:o + self.patch]
        x = (np.clip(bruto, self.lo, self.hi) - self.mu) / self.sd
        x[self.bin] = bruto[self.bin]                      # swir_estimado fica 0/1
        x = np.nan_to_num(x, nan=0.0)
        if self.aumentar:
            x = np.rot90(x, np.random.randint(4), axes=(1, 2))
            if np.random.rand() < 0.5:
                x = x[:, :, ::-1]
            ganho = 1 + 0.05 * np.random.randn(x.shape[0], 1, 1)       # outra cena
            desloc = 0.05 * np.random.randn(x.shape[0], 1, 1)
            ganho[self.bin], desloc[self.bin] = 1, 0
            x = x * ganho + desloc
        r = self.meta.iloc[j]
        return (torch.from_numpy(np.ascontiguousarray(x, dtype=np.float32)),
                torch.tensor(r.y, dtype=torch.float32), torch.tensor(r.peso, dtype=torch.float32))


# ═══════════════════════════════════════════════════════════════════════════
# MODELOS
# ═══════════════════════════════════════════════════════════════════════════
class SE(nn.Module):
    """Atenção de canal (squeeze-and-excitation)."""
    def __init__(self, c, r=8):
        super().__init__()
        self.f = nn.Sequential(nn.AdaptiveAvgPool2d(1), nn.Flatten(), nn.Linear(c, max(c // r, 4)),
                               nn.ReLU(inplace=True), nn.Linear(max(c // r, 4), c), nn.Sigmoid())

    def forward(self, x):
        return x * self.f(x)[:, :, None, None]


def bloco(ci, co, pool):
    camadas = [nn.Conv2d(ci, co, 3, padding=1, bias=False), nn.BatchNorm2d(co), nn.ReLU(inplace=True), SE(co)]
    return nn.Sequential(*camadas, nn.MaxPool2d(2)) if pool else nn.Sequential(*camadas)


class CNNPequena(nn.Module):
    def __init__(self, c_in):
        super().__init__()
        self.corpo = nn.Sequential(bloco(c_in, 32, False), bloco(32, 64, True), bloco(64, 64, False),
                                   bloco(64, 128, True), bloco(128, 128, False))
        self.cabeca = nn.Sequential(nn.AdaptiveAvgPool2d(1), nn.Flatten(), nn.Dropout(0.3), nn.Linear(128, 1))

    def forward(self, x):
        return self.cabeca(self.corpo(x)).squeeze(1)


class Backbone(nn.Module):
    def __init__(self, nome, c_in, patch, pretrained=True):
        super().__init__()
        import timm
        kw = {'img_size': patch} if nome.startswith('coatnet') else {}
        self.net = timm.create_model(nome, pretrained=pretrained, in_chans=c_in, num_classes=1,
                                     drop_rate=0.3, **kw)

    def forward(self, x):
        return self.net(x).squeeze(1)

    def congelar(self, sim):
        """Congela tudo menos a cabeça (classifier)."""
        cabeca = self.net.get_classifier()
        for p in self.net.parameters():
            p.requires_grad = not sim
        for p in cabeca.parameters():
            p.requires_grad = True


def criar_modelo(backbone, c_in, patch, pretrained=True):
    """pretrained=False na inferência: os pesos vêm do checkpoint."""
    return CNNPequena(c_in) if backbone == 'pequena' else Backbone(backbone, c_in, patch, pretrained)


# ═══════════════════════════════════════════════════════════════════════════
# TREINO E AVALIAÇÃO
# ═══════════════════════════════════════════════════════════════════════════
def fatores_classe(meta, idx):
    """Fatores que igualam o peso total de positivos (y ≥ 0,5) e negativos no treino."""
    pos = meta.y.values[idx] >= 0.5
    return 0.5 / max(pos.mean(), 1e-6), 0.5 / max((~pos).mean(), 1e-6)


def prever(modelo, loader, dispositivo):
    modelo.eval()
    saida = []
    with torch.no_grad(), torch.autocast(dispositivo.type, enabled=dispositivo.type == 'cuda'):
        for x, _, _ in loader:
            saida.append(torch.sigmoid(modelo(x.to(dispositivo, non_blocking=True))).float().cpu().numpy())
    return np.concatenate(saida)


def metricas(meta_v, p):
    m = meta_v.assign(p=p)
    rot = (m.y >= 0.5).astype(int)
    out: dict[str, float] = {'n': len(m)}
    if rot.nunique() == 2:
        out.update(auc_amostra=roc_auc_score(rot, m.p), ap=average_precision_score(rot, m.p))
    out['brier'] = float(np.mean((m.p - m.y) ** 2))
    oc = m[m.origem_amostra.str.startswith('ocorrencia')]
    neg = m[m.origem_amostra == 'negativo']
    if len(oc) and len(neg):
        ponto = oc.groupby('codigo_obj').p.median()
        out['auc_ponto'] = roc_auc_score(np.r_[np.ones(len(ponto)), np.zeros(len(neg))], np.r_[ponto, neg.p])
        out['n_pontos'] = len(ponto)
    return out


def treinar(X, meta, canais, idx_tr, idx_va, args, dispositivo, rotulo):
    """idx_va = None: treino final com número fixo de épocas, sem olhar validação."""
    est = estatisticas(X, canais, idx_tr, seed=args.seed)
    ds_tr = Patches(X, meta, idx_tr, canais, est, args.patch, aumentar=True)
    dl_tr = DataLoader(ds_tr, args.lote, shuffle=True, num_workers=args.workers, drop_last=True,
                       pin_memory=dispositivo.type == 'cuda')
    dl_va = None
    if idx_va is not None:
        ds_va = Patches(X, meta, idx_va, canais, est, args.patch, aumentar=False)
        dl_va = DataLoader(ds_va, args.lote * 2, shuffle=False, num_workers=args.workers)
    f_pos, f_neg = fatores_classe(meta, idx_tr)

    modelo = criar_modelo(args.backbone, len(canais), args.patch).to(dispositivo)
    pre = isinstance(modelo, Backbone)
    escalador = torch.amp.GradScaler(enabled=dispositivo.type == 'cuda')
    melhor, paciencia, hist, estado = -np.inf, 0, [], None

    def otimizador(lr):
        return torch.optim.AdamW([p for p in modelo.parameters() if p.requires_grad], lr=lr, weight_decay=1e-4)

    if pre:
        modelo.congelar(True)
    opt = otimizador(args.lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, args.epocas)
    for ep in range(args.epocas):
        if pre and ep == args.congelar:              # descongela tudo, lr 10× menor
            modelo.congelar(False)
            opt = otimizador(args.lr / 10)
            sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, args.epocas - ep)
        modelo.train()
        perdas = []
        t0 = time.time()
        for x, y, w in dl_tr:
            x, y, w = x.to(dispositivo), y.to(dispositivo), w.to(dispositivo)
            with torch.autocast(dispositivo.type, enabled=dispositivo.type == 'cuda'):
                logit = modelo(x)
                perda = F.binary_cross_entropy_with_logits(logit.float(), y, reduction='none')
            fator = torch.where(y >= 0.5, f_pos, f_neg)        # peso da amostra × equilíbrio de classes
            perda = (perda * w * fator).sum() / (w * fator).sum()
            opt.zero_grad(set_to_none=True)
            escalador.scale(perda).backward()
            escalador.step(opt)
            escalador.update()
            perdas.append(perda.item())
        sched.step()
        if dl_va is None:                                   # treino final: sem validação
            hist.append({'epoca': ep + 1, 'perda': float(np.mean(perdas))})
            print(f'  [{rotulo}] ép {ep + 1:3d}  perda {np.mean(perdas):.4f}  ({time.time() - t0:.0f} s)')
            continue
        p = prever(modelo, dl_va, dispositivo)
        m = metricas(meta.iloc[idx_va].reset_index(drop=True), p)
        alvo = np.nan_to_num(m.get('auc_ponto', m.get('auc_amostra', -m['brier'])), nan=-1e9)
        hist.append({'epoca': ep + 1, 'perda': float(np.mean(perdas)), **m})
        print(f'  [{rotulo}] ép {ep + 1:3d}  perda {np.mean(perdas):.4f}  '
              f'auc_amostra {m.get("auc_amostra", float("nan")):.3f}  '
              f'auc_ponto {m.get("auc_ponto", float("nan")):.3f}  ({time.time() - t0:.0f} s)')
        if alvo > melhor:
            melhor, paciencia = alvo, 0
            estado = {k: v.detach().cpu().clone() for k, v in modelo.state_dict().items()}
        else:
            paciencia += 1
            if paciencia >= args.paciencia:
                print(f'  [{rotulo}] parada antecipada na época {ep + 1}')
                break
    if dl_va is None:
        return modelo, est, pd.DataFrame(hist), None
    modelo.load_state_dict(estado)
    return modelo, est, pd.DataFrame(hist), prever(modelo, dl_va, dispositivo)


# ═══════════════════════════════════════════════════════════════════════════
def argumentos(lista=None):
    a = argparse.ArgumentParser(description=(__doc__ or '').split('\n')[1])
    a.add_argument('--modelo', choices=['A', 'B1'], default='A')
    a.add_argument('--backbone', default='pequena')
    a.add_argument('--patch', type=int, default=64)
    a.add_argument('--folds', default='0,1,2,3,4')
    a.add_argument('--teste', action='store_true',
                   help='treina em todos os folds (épocas fixas = --epocas) e avalia no teste')
    a.add_argument('--epocas', type=int, default=40)
    a.add_argument('--congelar', type=int, default=5, help='épocas só com a cabeça (backbones)')
    a.add_argument('--lote', type=int, default=64)
    a.add_argument('--lr', type=float, default=1e-3)
    a.add_argument('--paciencia', type=int, default=8)
    a.add_argument('--workers', type=int, default=4)
    a.add_argument('--seed', type=int, default=42)
    a.add_argument('--dados', default=str(PATCHES_DIR))
    a.add_argument('--saida', default=str(RUNS_DIR))
    return a.parse_args(lista)


def main(lista=None):
    args = argumentos(lista)
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    dispositivo = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    X, meta, canais = carregar(args.modelo, Path(args.dados))
    pasta = Path(args.saida) / f'{args.modelo}_{args.backbone}_{args.patch}px_{datetime.now():%Y%m%d_%H%M}'
    pasta.mkdir(parents=True, exist_ok=True)
    (pasta / 'config.json').write_text(json.dumps({**vars(args), 'canais': canais, 'n': len(meta),
                                                   'dispositivo': str(dispositivo)}, indent=2))
    print(f'{len(meta)} patches · {len(canais)} canais · {dispositivo} · saída {pasta}')
    print(pd.crosstab(meta.origem_amostra, meta.conjunto).to_string())

    treino = meta.conjunto == 'treino'
    if args.teste:
        # treina em todos os folds com --epocas fixas (use a época média das melhores do CV);
        # o teste não é olhado durante o treino
        idx_tr = np.where(treino)[0]
        idx_te = np.where(meta.conjunto == 'teste')[0]
        modelo, est, hist, _ = treinar(X, meta, canais, idx_tr, None, args, dispositivo, 'final')
        ds_te = Patches(X, meta, idx_te, canais, est, args.patch, aumentar=False)
        p = prever(modelo, DataLoader(ds_te, args.lote * 2, num_workers=args.workers), dispositivo)
        m = metricas(meta.iloc[idx_te].reset_index(drop=True), p)
        meta.iloc[idx_te].assign(p=p).to_csv(pasta / 'teste.csv', index=False)
        (pasta / 'metricas_teste.json').write_text(json.dumps(m, indent=2, default=float))
        torch.save({'estado': modelo.state_dict(), 'est': est, 'canais': canais, 'args': vars(args)},
                   pasta / 'final.pt')
        print(f'TESTE: {m}  (GTB: AUC por ponto {GTB_REF})')
        return

    linhas, oof = [], []
    for k in [int(f) for f in args.folds.split(',')]:
        idx_va = np.where(treino & (meta.fold == k))[0]
        idx_tr = np.where(treino & (meta.fold != k))[0]
        if not len(idx_va):
            continue
        print(f'\n═══ fold {k}: treino {len(idx_tr)}, validação {len(idx_va)}')
        modelo, est, hist, p = treinar(X, meta, canais, idx_tr, idx_va, args, dispositivo, f'fold{k}')
        hist.to_csv(pasta / f'historico_fold{k}.csv', index=False)
        torch.save({'estado': modelo.state_dict(), 'est': est, 'canais': canais, 'args': vars(args)},
                   pasta / f'fold{k}.pt')
        oof.append(meta.iloc[idx_va].assign(p=p))
        linhas.append({'fold': k, **metricas(meta.iloc[idx_va].reset_index(drop=True), p)})
    tab = pd.DataFrame(linhas)
    tab.to_csv(pasta / 'metricas_folds.csv', index=False)
    oof = pd.concat(oof, ignore_index=True)
    oof.to_csv(pasta / 'oof_validacao.csv', index=False)
    geral = metricas(oof, oof.p.values)
    print('\n' + tab.round(3).to_string(index=False))
    print(f'\nFORA DO FOLD (todos): {json.dumps(geral, default=float)}')
    print(f'Referência GTB: AUC por ponto {GTB_REF} (célula de fora, mosaico v1)')


if __name__ == '__main__':
    main()
