import matplotlib.pyplot as plt
import seaborn as sb
import numpy as np

from astropy.visualization import LogStretch, PercentileInterval
from code.crosstools import show_galaxy_rfgc, inellipse
import pandas as pd

from pathlib import Path
datapath = Path('data')
querypath = Path('queries')
name = f"rfgc_nearby_multiband2"
bands = 'grizy'
df = pd.read_hdf(datapath / (name + ".h5"))

for f in bands:
    df[f"{f}SerA"] = df[f"{f}SerRadius"]
    df[f"{f}SerB"] = df[f"{f}SerA"] * df[f"{f}SerAb"]
    df[f"{f}SerBa"] = 1/df[f"{f}SerAb"]

    df[f"{f}GalA"] = df[f"{f}GalMajor"]
    df[f"{f}GalB"] = df[f"{f}GalMinor"]
    df[f"{f}GalAb"] = df[f"{f}GalMinor"] / df[f"{f}GalMajor"]
    df[f"{f}GalBa"] = 1/df[f"{f}GalAb"]

df["AbO"] = df['bO'] / df['aO']
rfgc_only = df[
    [
        "RFGC",
        "PGC",
        "RAJ2000",
        "DEJ2000",
        "PA",
        "aO",
        "bO",
        "aE",
        "bE",
        "Btot",
        "AB",
        "MType",
        "Asym",
        "SB",
        "N",
    ]
].drop_duplicates()


ra, dec = df['raMean'], df['decMean']
ra_ref, dec_ref = df['RAJ2000'], df['DEJ2000']

ra = (ra - ra_ref)*np.cos(np.radians(dec_ref))
dec = dec - dec_ref

sigma_a = np.minimum(df['aO'], 10.)
sigma_b = np.minimum(df['bO'], 10.)

mask = inellipse((ra, dec), (0, 0), df['PA'], sigma_a, sigma_b)

# внутри контура галактики
matches_ref = df.loc[mask]

filt = 'g'

# есть galMajor/galMinor
wshape = matches_ref.dropna(subset=[f"{filt}GalIndex"])
# есть Серсик
fitted = wshape.dropna(subset=[f"{filt}SerRadius"])
# и там не чушь какая-то вписалась
sane = fitted[fitted[f"{filt}SerA"] > 0.051]
# и позиционный угол похож на тот что в RFGC
aligned = sane[
    (np.abs(sane[f"{filt}GalPhi"] - sane["PA"]) < 10)
  & (np.abs(sane[f"{filt}SerPhi"] - sane["PA"]) < 10)
]

# выбирается один наибольший по звездной величине источник
biggestsermag_aligned = (
    aligned.sort_values(["RFGC", f"{filt}SerMag"], ascending=[True, True])
    .groupby("RFGC")
    .head(1)
)

# настройка стиля
sb.set(rc={'figure.figsize': (4, 3)})
sb.set_style('whitegrid', {'grid.linestyle': ':'})
sb.set_palette("bright")

t = wshape
N = 2285

fig = plt.figure(figsize=(5, 5))
subset = t[t.RFGC == N]
profile_setup = {
    'sersic'    : True,
    'kron'      : False,
    'median'    : True,
    'petrosian' : False,
    'exp'       : False,
    'voculer'   : False
}
show_galaxy_rfgc(
    rfgc_only[rfgc_only.RFGC == N].iloc[0],
    filt,
    df = subset,
    fig = fig,
    image = True,
    zoom = 4,
    **profile_setup,
    transform = LogStretch(10) + PercentileInterval(99.5),
)
plt.show()
