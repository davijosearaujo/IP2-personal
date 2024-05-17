from Viscosity_models import *

num_modelos = 4

Tg = 500
m = 60
log_eta_inf = -4

liquidos = ['CaMgSi2O.xlsx', 'B2O3.xlsx', "NaAlSi3O8.xlsx"]

todas_as_equacoes = [VFT, AM, MYEGA, CW, F2, F4, F4, CLU, Bassler, FF, Demetriou_et_al]

equacoes = [VFT, AM, MYEGA, CW]
