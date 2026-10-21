### A. Strategy comparison (per-class promise 90%, share flagged on the random IID test set)
| Endpoint | Toxic rate | Random final | Most uncertain final | Uncertain+diverse final | Labels used by random | Labels needed by uncertain+diverse to match random's final |
|---|---|---|---|---|---|---|
| SR-ATAD5 | 3.7% | 0.870 | 0.797 | 0.855 | 1344 | 480 |
| SR-p53 | 6.2% | 0.837 | 0.700 | 0.747 | 1344 | 480 |
| NR-AhR | 11.7% | 0.387 | 0.439 | 0.428 | 1344 | not reached |
| SR-MMP | 15.8% | 0.438 | 0.476 | 0.452 | 1344 | not reached |

### Toxic compounds missed at the final round (mean over seeds)
| Endpoint | Standard, same chemistry | Per-class, same chemistry | Standard, new chemistry | Per-class, new chemistry |
|---|---|---|---|---|
| SR-ATAD5 | 99.6% | 3.6% | 100.0% | 1.3% |
| SR-p53 | 99.6% | 6.7% | 100.0% | 2.3% |
| NR-AhR | 79.6% | 9.0% | 92.2% | 7.5% |
| SR-MMP | 58.9% | 10.8% | 86.9% | 11.1% |

### B. Check-set size (miss rate on a random test set; promise = 10%)
| Toxic compounds in check set | Mean miss | 5th to 95th percentile | Theory mean | Theory 5th to 95th |
|---|---|---|---|---|
| 10 | 6.4% | 0.0% to 22.0% | 9.1% | 0.5% to 25.9% |
| 20 | 6.5% | 0.0% to 18.1% | 9.5% | 1.8% to 21.6% |
| 30 | 6.5% | 0.0% to 17.4% | 9.7% | 2.8% to 19.5% |
| 50 | 6.5% | 0.0% to 16.7% | 9.8% | 4.0% to 17.4% |
| 75 | 5.8% | 0.0% to 10.9% | 9.2% | 4.5% to 15.2% |

Random lab tests needed for 50 toxic compounds in the check set: SR-ATAD5 1,353, SR-p53 801, NR-AhR 426, SR-MMP 316

### C. Optional stopping (SR-p53, 30 toxic compounds in the check set, promise 10%)
- Mean miss at one pre-chosen round: 6.4%
- Mean miss when stopping the first time the workload target is met: 8.3% (above 15% in 19% of runs vs 2% for the pre-chosen round)
- Mean miss with 3 pre-chosen looks and the promise split 3 ways: 0.4%, but mean share flagged 99% vs 74% for the naive rule (the price of the safer rule)

### D. New chemistry
| Endpoint | Alarm detects | False alarm | Miss, same chem. | Miss, new chem. | Flagged, same chem. | Flagged, new chem. | Flagged after restart | Miss after restart | Restart lab tests |
|---|---|---|---|---|---|---|---|---|---|
| SR-ATAD5 | 100% | 11% | 5.8% | 4.3% | 86% | 94% | 91% | 8.4% | 722 |
| SR-p53 | 100% | 1% | 11.6% | 8.6% | 69% | 90% | 86% | 9.4% | 690 |
| NR-AhR | 100% | 6% | 5.8% | 4.6% | 49% | 81% | 70% | 10.0% | 668 |
| SR-MMP | 100% | 5% | 6.9% | 4.8% | 59% | 87% | 78% | 8.1% | 591 |

### A2. Plateau stop rule (stop when workload improved by under 1 point over the last 2 rounds; needs 4+ rounds and workload under 98%)
| Endpoint | Strategy | Mean labels at stop | Mean share flagged at stop | Share flagged if run to the end |
|---|---|---|---|---|
| SR-ATAD5 | Random | 1306 | 86.7% | 87.0% |
| SR-ATAD5 | Most uncertain | 998 | 85.3% | 79.7% |
| SR-ATAD5 | Uncertain + diverse | 979 | 80.2% | 85.5% |
| SR-p53 | Random | 1018 | 80.1% | 83.7% |
| SR-p53 | Most uncertain | 768 | 76.7% | 70.0% |
| SR-p53 | Uncertain + diverse | 768 | 76.3% | 74.7% |
| NR-AhR | Random | 730 | 50.0% | 38.7% |
| NR-AhR | Most uncertain | 749 | 53.4% | 43.9% |
| NR-AhR | Uncertain + diverse | 749 | 52.8% | 42.8% |
| SR-MMP | Random | 730 | 61.2% | 43.8% |
| SR-MMP | Most uncertain | 1037 | 56.0% | 47.6% |
| SR-MMP | Uncertain + diverse | 845 | 58.6% | 45.2% |