# MeteoD S2S ECDS → WebP → Cloudflare R2 (v10 – soglie ERA5 su griglia)

**Senza Leonardo.** Ogni giorno GitHub Actions scarica da ECMWF ECDS la previsione S2S della **data UTC di due giorni prima**. Interpola i membri dell'ensemble sulla griglia **0,4°** degli indici `assets/province_index_mediterraneo_0.4.json`; calcola i quintili usando i percentili **ERA5 P20/P40/P60/P80 per ogni pixel** da `assets/era5_week_percentiles_2006_2025_conservative_04_smooth.nc` e genera i plot MeteoD in **WebP**. Alla fine li carica su Cloudflare R2. Non crea NetCDF, JSON di risultati, né altri output da pubblicare. I GRIB intermedi sono in una cartella temporanea cancellata automaticamente.


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
**Nessuna correzione viene applicata alle precipitazioni.** Le soglie cittadine storiche nel JSON non sono modificate, ma non vengono più
utilizzate per calcolare le probabilità.

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

## Soglie dei quintili: ERA5 2006–2025, griglia 0,4° (v10)

L'asset statico `assets/era5_week_percentiles_2006_2025_conservative_04_smooth.nc`
contiene i **quattro percentili** P20, P40, P60, P80 per `t2m` (°C) e
`tp` (mm per settimana), ciascuno con 52 mappe geografiche. Questo è
**l'unico input delle soglie** per classificare gli ensemble, sia nelle
mappe sia negli istogrammi per le 151 località.

- All'avvio, i campi ERA5 vengono letti e controllati (variabili, coordinate,
  unità, NaN, ordine P20 ≤ P40 ≤ P60 ≤ P80). Un file mancante o incompleto
  fa interrompere l'esecuzione **prima del download ECDS**.
- I dati sono su **centri di cella** (58×135, incluso l'ultimo mezzo intervallo
  meridionale), mentre i campi previsionali usano 58×136 nodi. Le mappe dei
  percentili sono interpolate **bilinearmente** sui nodi della griglia
  previsionale; solo per i punti appena oltre il primo/ultimo centro di cella
  si utilizza il valore del centro più vicino (nessuna estrapolazione distante).
- La climatologia ha 52 settimane di 7 giorni fissi a partire dal 1° gennaio.
  La previsione è lunedì–domenica: per ciascuna settimana previsionale si
  applica una media delle **sette mappe di soglia associate ai sette giorni**.
  Questo è coerente con l'allineamento temporale del bias, ma resta
  un'approssimazione: i percentili della settimana lunedì–domenica esatta
  non si possono ricostruire rigorosamente dai soli percentili di settimane
  fisse.
- Per ciascun pixel, i 50 membri (temperatura già corretta per il bias,
  precipitazione non corretta) vengono confrontati con le quattro soglie.
  Quindi `Q1: <=P20`, `Q2: (P20,P40]`, `Q3: (P40,P60]`,
  `Q4: (P60,P80]`, `Q5: >P80`. La probabilità è la percentuale
  di membri che ricadono nel quintile.
- **Le stesse soglie grigliate** vengono lette nel pixel della città per
  calcolare le percentuali degli istogrammi. In particolare, i numeri
  in `assets/soglie_settimanali_med.json.gz` NON influenzano più la
  classificazione (neppure sulle mappe). Il vecchio JSON è conservato
  **esclusivamente per stampare P0/P20/P40/P60/P80/P100 ai bordi delle
  barre degli istogrammi**, come esplicitamente richiesto. Per questo tali
  etichette possono differire dalle soglie grigliate che determinano le
  percentuali visualizzate.

La palette per le mappe (Q1+Q2 / Q3 / Q4+Q5, peso doppio al centrale e
contorni Q1/Q5 ≥80%) non è stata modificata. Il download congiunto di
`2_m_temperature` e `total_precipitation`, la correzione del bias di tutti
i membri della sola temperatura, gli istogrammi e la pubblicazione R2
rimangono invariati.

## Directory di Cloudflare confermata dal vecchio codice

Il codice Leonardo aveva `R2_BUCKET_NAME = "settimanali"` e chiamava `upload_file(file_path, "settimanali", filename)`. Quindi il percorso è la **radice del bucket `settimanali`**, non `settimanali/PLOT/` né altre sottocartelle. La pipeline replica questa disposizione: `R2_BUCKET=settimanali` e `R2_PREFIX=''`. È possibile impostare un prefisso diverso soltanto se il sito usa realmente un'altra directory.

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

Gli unici file di output locali sono sotto `output/20261022/mappe/` e `output/20261022/citta/`, tutti `.webp`. Non esistono JSON di output, NetCDF o matrici salvate. Gli asset NetCDF e JSON inclusi sono **input statici**: il NetCDF dei percentili determina le probabilità, il NetCDF del bias corregge la temperatura, il JSON delle soglie cittadine stampa soltanto i valori ai confini dei quintili negli istogrammi.

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
- **Risolto in v10:** sia le mappe sia gli istogrammi usano esclusivamente soglie ERA5 calcolate pixel per pixel; non si interpolano più le soglie dalle città. Le etichette numeriche sui bordi degli istogrammi provengono invece dal vecchio JSON e potrebbero non coincidere esattamente con le soglie effettive.
- Le aree colorate continuano a usare l'approssimazione grafica originale `40/20/40` riscalata, **non veri tercili**.
- L'accesso ECMWF reale, la decodifica dei GRIB scaricati, le mappe Cartopy e l'upload reale R2 necessitano della prima esecuzione completa su GitHub. I test offline verificano calcoli e WebP e simulano l'upload senza trasferire dati.
- L'upload R2 non è una transazione atomica: se fallisce dopo alcuni oggetti, il workflow termina in errore e alcuni file potrebbero risultare già aggiornati; basta rilanciare la run con la stessa data.
- Per ora le soglie del file coprono il 2026–2030. Quando non sono disponibili per una data, la procedura si ferma senza pubblicare risultati incompleti.

## Dettagli ECMWF

Le richieste ora usano **esattamente i nomi dei campi del form ECDS**: `origin`, `year`, `month`, `day`, `time`, `level_type`, `variable`, `forecast_type`, `leadtime_hour`, `data_format`, `area`. Non si usa più la richiesta MARS manuale con `class`, `param` e `step`.

Si effettua **una sola richiesta ECDS** al dataset `s2s-forecasts` con
`variable=["2_m_temperature", "total_precipitation"]` e la selezione congiunta dei 92 valori `leadtime_hour`:

- temperatura: 46 intervalli giornalieri `0_24`, `24_48`, ..., `1080_1104`;
- precipitazione: 46 accumuli progressivi alle ore `24`, `48`, ..., `1104`.

ECDS prepara un solo job e restituisce un unico download. La pipeline separa
**localmente** i messaggi GRIB usando i metadati fisici (`shortName`, categoria,
stepType), non la posizione nel file. Accetta sia uno stream GRIB sia un archivio ZIP
contenente GRIB; li passa separatamente a `cfgrib`, che gestisce le diverse
coordinate temporali di temperatura e precipitazione. I GRIB intermedi restano
nel `TemporaryDirectory` e sono eliminati dopo l'elaborazione. La precipitazione
cumulata viene ancora differenziata per ottenere gli apporti giornalieri (mm).

L'API reale ECDS non è verificabile nei test offline: se il backend non accetta
questa selezione mista di lead time, la pipeline **fallisce esplicitamente** e
non scarica i due parametri con due richieste di ripiego. Anche quando il download
riesce, vengono controllati i 46 intervalli giornalieri e la coerenza dei membri.
 Scaricare anche i 6-hourly non serve per questa applicazione. La richiesta usa l'area ampliata `[54, -14, 28, 44]` esclusivamente per mantenere i punti necessari all'interpolazione bilineare: l'output finale resta sulla griglia 0,4° `[52, -12, 29, 42]`. Le longitudini GRIB 0–360° vengono ricondotte a −180–180° prima di interpolare.

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


## v9 — Download unico ECMWF ECDS (temperatura + precipitazione)

Il download è ora un solo job per data (`ecmwf_s2s_t2m_tp_YYYYMMDD.grib`),
contenente entrambe le variabili e l'unione dei loro lead time. Dopo il download
lo script genera **in locale** due GRIB temporanei (`ecmwf_s2s_t2m_*.grib` e
`ecmwf_s2s_tp_*.grib`), separando i messaggi con ecCodes. Questo permette di
riutilizzare gli stessi metodi di aggregazione: temperatura media dei 7 giorni,
precipitazione accumulata dei 7 giorni. La correzione di bias continua ad essere
applicata pixel per pixel a ogni membro **solo alla temperatura**.

`python pipeline.py --dry-run` stampa la singola richiesta ECDS, senza scaricare.
`--grib-dir` accetta sia il nuovo file congiunto sia la coppia di file precedenti.

Il supporto della combinazione di intervalli e lead time cumulativi in un'unica
richiesta deve essere verificato con la prima esecuzione reale sul servizio ECDS.
