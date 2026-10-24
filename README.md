# MeteoD - Settimanali ECMWF S2S, v2 (senza Leonardo)

## Cosa fa
1. Scarica **temperature medie giornaliere** e **precipitazioni cumulative** dell'**ensemble perturbato ECMWF S2S** dal portale **ECDS** (dataset `s2s-forecasts`).
2. Ottiene gli accumuli quotidiani di pioggia per differenza delle precipitazioni cumulative (S2S `tp`); la temperatura viene convertita da K a °C.
3. Converte in **5 settimane lunedì-domenica complete**, partendo dal primo lunedì non precedente alla run. Temperatura: media sui sette giorni; precipitazione: totale sui sette giorni.
4. Esegue **interpolazione bilineare a 0,4°** prima del confronto con le soglie. Griglia di **58 x 136** punti, latitudini **52,0; 51,6; …; 29,2** (decrescenti) e longitudini **−12; −11,6; …; 42** (crescenti), mantenendo invariati gli indici del file `province_index_mediterraneo_0.4.json`.
5. Calcola i cinque quintili per ciascuna località usando i **valori originali** di `soglie_settimanali_med.pkl`, qui esportato nel file sicuro `assets/soglie_settimanali_med.json.gz` (151 località, 257 settimane da febbraio 2026 a dicembre 2030).
6. Produce gli **stessi tipi di plot** del vecchio script Leonardo, con la parte grafica preservata in `plotting.py`: istogramma temperatura/precipitazione per località (prime quattro settimane); mappa per ogni settimana e variabile (5×2 WebP).
7. Scrive il JSON delle località nel vecchio formato `T_1..T_5`, `P_1..P_5` e un NetCDF con i valori dell'ensemble settimanale interpolato.

### Importante: limitazione climatologica per le MAPPE
Il file di soglie fornito contiene valori **solo per località**, non una climatologia ERA5 di griglia. Le mappe v1 utilizzano dunque la **interpolazione spaziale lineare delle soglie tra le 151 località**, con assegnazione dal punto più vicino fuori dall'inviluppo delle località. Località con identico indice sono mediate **solo per costruire le mappe**. Le località e i relativi istogrammi mantengono **ognuna le proprie soglie originali**, senza mediare o interpolare soglie tra città.

Le mappe v1 sono quindi una **APPROSSIMAZIONE**: mantengono stile e filename Leonardo, ma non possono coincidere scientificamente con le mappe originali. La fase successiva consigliata sarà una **climatologia ERA5 grigliata** a 0,4°, in modo da poter calcolare quintili senza interpolare soglie tra città.

Per replicare l'aspetto delle mappe viene conservato il vecchio trucco grafico dei tre gruppi `Q1+Q2`, `Q3`, `Q4+Q5` con riscalatura del gruppo centrale: **non sono tercili meteorologici veri**.

Non vengono effettuati downscaling statistico, correzione bias o calibrazione su reforecast: è la v1 sperimentale, da confrontare con i plot Leonardo.

## Primo setup su GitHub

1. Crea un nuovo repository (preferibilmente **privato**) e carica **il contenuto** di questa cartella, compresa `.github/workflows/subseasonal.yml` e `assets/`. Non caricare il vecchio script con chiavi R2 dentro.
2. Accedi a https://ecds.ecmwf.int/datasets/s2s-forecasts?tab=download e **accetta le condizioni/licenze** del dataset.
3. Vai su https://ecds.ecmwf.int/how-to-api e recupera il tuo token personale ECMWF.
4. Su GitHub vai in **Settings → Secrets and variables → Actions → New repository secret** e crea `ECDS_API_KEY` con quel token. **Non** inserirlo negli script.
5. Vai in **Actions → MeteoD ECMWF S2S Settimanali → Run workflow**. Per la prima prova lascia il campo `run_date` **vuoto**. Verrà selezionata automaticamente la data UTC di **due giorni prima**, indipendentemente dal giorno della settimana. Per un test storico puoi inserire ad esempio `2026-10-22`.
6. Verifica i log. Nel riepilogo della run scarica l'artifact `meteod-s2s-output-N`, con WebP, JSON e NetCDF.

È presente una pianificazione **tutti i giorni alle 12:15 UTC**, che richiede la run ECMWF S2S delle 00 UTC di due giorni prima (circa 60 ore dopo l'inizializzazione). ECMWF S2S è giornaliero dall'aggiornamento di novembre 2024. Nessun upload su Cloudflare R2 per questa **prima versione di validazione**.

### Esecuzione locale

Richiede Python 3.11/3.12, librerie in `requirements.txt` ed ecCodes per cfgrib. In Linux, esegui prima `sudo apt-get install libeccodes0 libeccodes-dev`.

```bash
pip install -r requirements.txt
export ECDS_API_KEY='IL_TUO_TOKEN'
python pipeline.py --dry-run
python pipeline.py --output output
# Per una emissione storica specifica:
python pipeline.py --date 2026-10-22 --output output
```

Per ripetere senza riscaricare GRIB già presenti, la pipeline salta i file non vuoti. Per usare file GRIB propri: `python pipeline.py --date 2026-10-22 --grib-dir MIA_CARTELLA` con nomi `ecmwf_s2s_t2m_20261022.grib` e `ecmwf_s2s_tp_20261022.grib`.

### Output attesi

```text
output/20261022/
├── grib/
│   ├── ecmwf_s2s_t2m_20261022.grib
│   └── ecmwf_s2s_tp_20261022.grib
├── ensemble_settimanale_04.nc
├── province_quintili_2026_10_22.json
├── METODOLOGIA.txt
├── PLOT/
│   ├── subseasonal_20261026_2tp.webp
│   ├── subseasonal_20261026_tpp.webp
│   ├── ... (cinque settimane, 10 mappe)
└── PLOT_istogrammi/
    ├── Bologna_20261026.webp
    ├── Milano_20261026.webp
    └── ... (151 località con 4 settimane per grafico)
```

### Test

```bash
python tests/test_offline.py
```

È un test offline con 12 membri sintetici, nessuna chiamata ECMWF. Verifica differenza delle precipitazioni cumulative, griglia, interpolazione, cinque quintili, tutte le località, struttura delle mappe, NetCDF e un istogramma WebP. Il download reale ECDS, i GRIB scaricati e le mappe Cartopy complete **devono ancora essere verificati alla prima esecuzione Actions**.

### Modifiche v2: data giornaliera e errore MARS `Ambiguous`

* Se non viene passato `--date`, la pipeline usa la data **UTC di due giorni prima**; dal 12 novembre 2024 ECMWF S2S è giornaliero, non solo lunedì e giovedì.
* La richiesta ECDS viene costruita nel formato **MARS-like** con codici numerici: `class=27` (S2S), `param=167` (temperatura 2 m), `param=228` (precipitazione totale). Questo evita di inviare i nomi descrittivi ambigui che sono comparsi nel log `MarsRuntimeError`.
* Dopo il download tutte le altre funzioni della v1 (interpolazione 0,4°, province, soglie, mappe e istogrammi) rimangono inalterate.
* **Importante**: i codici sono quelli ufficiali ECMWF, ma non è stato possibile provare dal vivo che il backend ECDS accetti `class=27` nella sintassi MARS-like. Se si riceve un errore diverso (ad es. class non valida), il log della prima esecuzione sarà necessario per adeguare la richiesta. Non attribuire eventuali errori di sintassi alla mancanza dei dati.
* Il ritardo di 48 ore è un vincolo minimo, **non una garanzia di pubblicazione** a un preciso minuto. Se il file richiesto non è stato ancora archiviato, il job fallisce e non utilizza silenziosamente una data precedente.

### Errori comuni

- `403` o problemi di autorizzazione: accetta la licenza ECDS con l'account associato al token.
- `MarsRuntimeError`/`Ambiguous`: nella v2 i nomi ambigui sono stati sostituiti da codici numerici. Se persiste, confrontare i log con la richiesta MARS-like generata dal portale ECDS.
- `Step incompleti`: il forecast ottenuto non contiene gli intervalli 24h previsti; non proseguire con aggregazioni parziali.
- `Precipitazione cumulata in diminuzione`: fermarsi e controllare i metadati del GRIB, nessuna conversione automatica cieca.
- `Soglie non disponibili`: il file contiene solo settimane **2026–2030**. Nuove date richiedono l'aggiornamento della climatologia di riferimento.
- I pacchetti Cartopy possono scaricare file Natural Earth (coste e confini) alla prima creazione delle mappe, richiedendo accesso Internet.

## Fonti tecniche
- https://ecds.ecmwf.int/how-to-api
- https://ecds.ecmwf.int/datasets/s2s-forecasts?tab=overview
- https://confluence.ecmwf.int/spaces/S2S/pages/26903293/Parameters
- https://confluence.ecmwf.int/spaces/S2S/pages/27399606/Total+precipitation
