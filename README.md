# MeteoD S2S ECDS → WebP → Cloudflare R2 (v6)

**Senza Leonardo.** Ogni giorno GitHub Actions scarica da ECMWF ECDS la previsione S2S della **data UTC di due giorni prima**. Interpola i membri dell'ensemble sulla griglia **0,4°** degli indici `assets/province_index_mediterraneo_0.4.json`; calcola i quintili rispetto alle soglie di `assets/soglie_settimanali_med.json.gz` e genera i plot MeteoD in **WebP**. Alla fine li carica su Cloudflare R2. Non crea NetCDF, JSON di risultati, né altri output da pubblicare. I GRIB intermedi sono in una cartella temporanea cancellata automaticamente.


## Correzione bias settimanale della temperatura (nuova versione)

L'asset **`assets/bias_s2s_minus_era5_2015_2020_04_smooth.nc`** contiene il
bias climatologico della temperatura media a 2 m nel periodo 2015–2020:
`bias = S2S t2m - ERA5 t2m` (52 settimane, °C). È un file di **input statico**
incluso nel repository: non viene caricato su Cloudflare e non deve essere
scaricato a ogni esecuzione. Per sostituirlo con una climatologia aggiornata,
sostituire il file in `assets` senza cambiare codice.

Dopo l'interpolazione dei GRIB sulla griglia operativa 0,4°, ma **prima** di
calcolare probabilità/istogrammi, ogni membro dell'ensemble è corretto in ogni
pixel: `T_corretta = T_S2S - bias`. La stessa correzione additiva è applicata
a tutti i membri di un pixel; preserva quindi la dispersione ensemble.
**Nessuna correzione viene applicata alle precipitazioni.** Le soglie cittadine
esistenti non sono modificate.

Il bias usa settimane climatologiche fisse dall'1 gennaio (1–7 gennaio,
8–14 gennaio, ..., settimana 52). Le settimane previsionali sono lunedì–domenica
e possono attraversare due settimane climatologiche o il nuovo anno: il bias
per la previsione è la media **dei 7 bias giornalieri corrispondenti** alla
finestra lunedì–domenica (bias giornaliero inteso come mappa della relativa
settimana climatologica). Il 29 febbraio viene assegnato alla settimana del
28 febbraio e il 31 dicembre alla settimana 52, perché la climatologia
contiene soltanto 52 × 7 giorni.

Il NetCDF fornito ha coordinate **centro-cella** di 58×135 punti, mentre la
griliate di previsione originale usa 58×136 **nodi**. Il bias viene quindi
riallineato una sola volta con interpolazione bilineare. Nei dati allegati sono
presenti 3 righe di bias mancanti a nord e 3 a sud (in ogni settimana): tali
zone vengono riempite con il pixel disponibile più vicino prima del
riallineamento. Anche i nodi appena esterni ai centri delle celle sono limitati
al centro valido più vicino; la correzione in quelle fasce costituisce quindi
una **estrapolazione approssimata**, non un bias osservato localmente.
Il log segnala quanti valori sono stati estesi. Se il file manca o le
coordinate/unità/definizione sono incompatibili, la pipeline si ferma prima
del download.

Le mappe di probabilità e i 151 istogrammi delle località useranno le
**temperature corrette**. I prodotti pubblicati rimangono solo WebP e le
precipitazioni non cambiano. Non vengono prodotti NetCDF di output.

## Directory di Cloudflare confermata dal vecchio codice

Il codice Leonardo aveva `R2_BUCKET_NAME = "settimanali"` e chiamava `upload_file(file_path, "settimanali", filename)`. Quindi il percorso è la **radice del bucket `settimanali`**, non `settimanali/PLOT/` né altre sottocartelle. La v6 replica questa disposizione: `R2_BUCKET=settimanali` e `R2_PREFIX=''`. È possibile impostare un prefisso diverso soltanto se il sito usa realmente un'altra directory.

Esempi di chiavi R2:

```text
settimanali/subseasonal_20261026_2tp.webp
settimanali/subseasonal_20261026_tpp.webp
settimanali/Bologna_20261026.webp
settimanali/Milano_20261026.webp
```

**Formato:** questa versione torna al formato `.webp` del vecchio flusso Leonardo. Le mappe e i grafici locali sono salvati con `dpi=150`, `bbox_inches="tight"` e `quality=70`. L'upload usa `Content-Type: image/webp`. Gli eventuali vecchi file `.png` su R2 non vengono cancellati automaticamente.

## Configurazione GitHub

1. Caricare nel repository il **contenuto** della directory (compresa `.github/workflows/`). Non caricare vecchi sorgenti che contengono chiavi cloud in chiaro.
2. In **GitHub → Settings → Secrets and variables → Actions → New repository secret**, creare:
   - `ECDS_API_KEY` — token personale ECMWF ECDS;
   - `R2_ACCOUNT_ID` — account ID Cloudflare;
   - `R2_ACCESS_KEY_ID` — Access Key di una credenziale R2 con autorizzazione di scrittura sul bucket `settimanali`;
   - `R2_SECRET_ACCESS_KEY` — Secret Key associata.
3. Rigenerare/revocare le vecchie chiavi R2 che apparivano nello script Leonardo e **non riutilizzarle**.
4. Accettare la licenza del dataset su https://ecds.ecmwf.int/datasets/s2s-forecasts?tab=download.
5. In **Actions → MeteoD ECMWF S2S → Cloudflare R2 → Run workflow**, lasciare `run_date` vuoto. Il cron automatico gira tutti i giorni alle **12:15 UTC**.

## Output

Per una data ECMWF come `2026-10-22` genera **161 WebP**:

- 10 mappe: 5 settimane × temperatura (`2tp`) e precipitazione (`tpp`).
- 151 istogrammi per le località del file `province_index_mediterraneo_0.4.json`, ciascuno con le prime quattro settimane.

Gli unici file di output locali sono sotto `output/20261022/mappe/` e `output/20261022/citta/`, tutti `.webp`. Non esistono JSON di output, NetCDF o matrici salvate. I riferimenti `assets/*.json*` sono **input statici indispensabili** per costruire gli istogrammi e i campi di probabilità.

Il caricamento sul bucket mantiene i soli nomi dei file, senza la struttura `mappe/` e `citta/`, come faceva Leonardo. L'header HTTP `Content-Type` è `image/webp` e il TTL cache è 300 secondi, utile perché gli oggetti per una stessa settimana possono venire aggiornati il giorno dopo. Le mappe future vengono sovrascritte dalla previsione più recente se il loro nome coincide; gli oggetti più vecchi non vengono cancellati.

## Esecuzione locale senza upload

```bash
python -m pip install -r requirements.txt
python tests/test_offline.py
python pipeline.py --dry-run
python pipeline.py --date 2026-10-22 --no-upload --output output
```

L'esecuzione normale `python pipeline.py --output output` esegue anche l'upload su Cloudflare e richiede i secrets come variabili ambiente. `--grib-dir DIR` può leggere GRIB già scaricati, senza eliminarli. `--no-maps` genera soltanto le immagini delle località.

## Limitazioni della versione semplificata

- **Aggiornato:** ora la temperatura viene corretta mediante il NetCDF del bias prima di mappe e istogrammi, mentre le precipitazioni mantengono la procedura originale.
- Le soglie fornite sono per *località*. Le mappe utilizzano soglie interpolate spazialmente e sono quindi **approssimazioni**, non quintili ottenuti da una climatologia ERA5 grigliata.
- Le aree colorate continuano a usare l'approssimazione grafica originale `40/20/40` riscalata, **non veri tercili**.
- L'accesso ECMWF reale, la decodifica dei GRIB scaricati, le mappe Cartopy e l'upload reale R2 necessitano della prima esecuzione completa su GitHub. I test offline verificano calcoli e WebP e simulano l'upload senza trasferire dati.
- L'upload R2 non è una transazione atomica: se fallisce dopo alcuni oggetti, il workflow termina in errore e alcuni file potrebbero risultare già aggiornati; basta rilanciare la run con la stessa data.
- Per ora le soglie del file coprono il 2026–2030. Quando non sono disponibili per una data, la procedura si ferma senza pubblicare risultati incompleti.

## Dettagli ECMWF

Le richieste ora usano **esattamente i nomi dei campi del form ECDS**: `origin`, `year`, `month`, `day`, `time`, `level_type`, `variable`, `forecast_type`, `leadtime_hour`, `data_format`, `area`. Non si usa più la richiesta MARS manuale con `class`, `param` e `step`.

Si effettuano **due richieste** al dataset `s2s-forecasts` per evitare che i lead times delle due variabili siano mescolati:

- `variable=["2_m_temperature"]`: i 46 intervalli giornalieri `0_24`, `24_48`, ..., `1080_1104`;
- `variable=["total_precipitation"]`: i 46 accumuli progressivi alle ore `24`, `48`, ..., `1104`. Vengono differenziati per ottenere il totale di ogni singolo giorno (mm). Scaricare anche i 6-hourly non serve per questa applicazione. La richiesta usa l'area ampliata `[54, -14, 28, 44]` esclusivamente per mantenere i punti necessari all'interpolazione bilineare: l'output finale resta sulla griglia 0,4° `[52, -12, 29, 42]`. Le longitudini GRIB 0–360° vengono ricondotte a −180–180° prima di interpolare.

L'accesso reale a ECDS non è stato verificato offline: **non è possibile garantire che l'errore interno `MarsRuntimeError` non si ripresenti** anche passando dalla sintassi ECDS ufficiale.

Riferimenti:
- https://ecds.ecmwf.int/datasets/s2s-forecasts?tab=download
- https://confluence.ecmwf.int/spaces/S2S/pages/27399606/Total+precipitation
- https://confluence.ecmwf.int/spaces/DAC/pages/540541961/Examples+of+the+current+and+future+retrieval+requests


Riferimenti: https://ecds.ecmwf.int/how-to-api e https://confluence.ecmwf.int/spaces/DAC/pages/540541961/Examples+of+the+current+and+future+retrieval+requests.


## v5 — Correzione GitHub exit 134 dopo upload

Se nel log GitHub compare `DONE: 161 WebP caricati su R2` seguito da
`double free or corruption (!prev)` / `Aborted (core dumped)` (exit 134),
si tratta di un arresto del processo Python, **non necessariamente di un
fallimento del trasferimento**. Il crash viene da codice nativo; non è
possibile attribuirlo con certezza a pyproj, Cartopy, ecCodes o altre librerie
senza un backtrace.

Per mitigare questo problema, il workflow **v5**:

1. scarica e genera i WebP come in v4;
2. carica tutti gli oggetti in Cloudflare R2;
3. esegue `head_object()` per **ogni** WebP e confronta la dimensione remota con quella locale (impostazione `R2_VERIFY_OBJECTS=1`);
4. **soltanto se tutto è terminato senza eccezioni**, salta la fase di
   finalizzazione dell'interprete Python con `os._exit(0)` nel comando GitHub,
   evitando il crash osservato al termine del processo.

La chiamata `os._exit(0)` è una **mitigazione mirata al crash di chiusura**,
non una correzione della causa nativa: non esegue i distruttori Python
restanti. È usata nel solo workflow dopo il completamento di tutte le
operazioni; il codice `pipeline.py` eseguito normalmente mantiene il normale
ciclo di vita e può ancora mostrare il crash durante debug.

`PYTHONFAULTHANDLER=1` aiuta a ottenere informazioni se un futuro crash si
verifica *durante* l'elaborazione. Gli avvisi `pyproj unable to set PROJ
database path` e `Cartopy Downloading Natural Earth` vanno distinti da un
crash; l'origine dell'errore nativo resta da indagare, per esempio fissando
versioni note delle dipendenze o isolando la generazione mappe dalla lettura
GRIB in processi diversi.


## v6 — Output WebP (24 ottobre 2026)

- Tutti i 161 grafici sono `.webp`: 10 mappe e 151 istogrammi.
- Entrambe le funzioni Matplotlib usano `format="webp"`, `dpi=150`, `bbox_inches="tight"`, `pil_kwargs={"quality": 70}`.
- L'upload cerca esclusivamente `.webp`, controlla il numero atteso di file, imposta `ContentType="image/webp"` e verifica le dimensioni remote con `head_object()` se `R2_VERIFY_OBJECTS=1`.
- La destinazione resta **la radice** del bucket R2 `settimanali`. I PNG precedenti rimangono nel bucket fino all'eventuale cancellazione manuale.
- Rimane attiva la mitigazione v5 del crash nativo in chiusura Python, usata solo nel workflow GitHub dopo upload e verifica completati.


## Aggiornamento v7

Le mappe WebP delle settimanali ora seguono il blocco di plotting del flusso originario Leonardo, inclusi layout, legenda, hatching degli eventi estremi e naming dei file `subseasonal_YYYYMMDD_{2tp|tpp}.webp`.
