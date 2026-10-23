# CSDA V3 - istruzioni rapide

Il pacchetto riproduce i risultati salvati del manoscritto CSDA senza rifittare automaticamente gli stimatori. La Table S4 usa `results/reference/fixed32/FIXED_all_results.csv`: 32 dataset, 6 metodi, 192 fit; REG-HD-TFKM usa HIER4 e REG-HD-TRKM usa RV--CV.

Da MATLAB:
```matlab
RUN_REPLICATION("CHECK");
RUN_REPLICATION("REPORT");
RUN_REPLICATION("SELFTEST");
```
I microdati bancari confidenziali non sono distribuiti; il report pubblico usa solo aggregati non identificativi.
