# When can the model stop asking the lab?

A retrospective simulation of an active learning loop for toxicity screening on the public Tox21 data.
It tests one idea: let the model choose which compounds go to the lab next, and stop only when a
**per-class promise** about the rare toxic compounds can be backed by a **randomly drawn check set**.

## The question
Each lab test costs money. A model can rank compounds, but when is it safe to stop testing?
This project asks four concrete things:
1. Does choosing batches smartly (instead of randomly) cut the number of lab tests?
2. Does a per-class (Mondrian) conformal promise keep the miss rate on toxic compounds, where a standard overall promise does not?
3. How many lab tests must go into the check set, and what does stopping "the first time it looks good" do to the promise?
4. Can a simple similarity alarm tell when new chemistry arrives, and what does restarting buy?

## What it does
- Data: Tox21 (public), four endpoints with 4%, 6%, 12% and 16% toxic compounds. Morgan fingerprints, random forest.
- Splits: a random check set (calibration), a random test set, and a "new chemistry" test set made of rare scaffolds the loop never sees.
- Loop: start with 192 random labels, then 12 rounds of 96 compounds (plate size). Strategies: random, most uncertain, uncertain + diverse.
- Promise: "we flag at least 90% of toxic compounds" (miss at most 10%), using split conformal prediction calculated for the toxic class alone.
- The check set is chosen at random and kept apart from the actively selected training set. Its lab cost is counted.

## Results (5 seeds for A; see results/summary.md for every number)
1. **A standard promise hides the miss rate on the rare class.** At the final round, a standard (overall) conformal promise missed 59% to 100% of toxic compounds. The per-class promise missed 4%, 7%, 9% and 11% (target: at most 10%; SR-MMP was slightly over). This matches a July 2026 arXiv preprint (not peer reviewed).
2. **Smart batch selection helped on the two rare-toxic endpoints, not on the two common ones.** On SR-ATAD5 and SR-p53, uncertainty-based selection ended with fewer compounds flagged than random (80% vs 87%, and 70% vs 84%). On NR-AhR and SR-MMP, random was as good or better. Adding diversity did not beat plain uncertainty.
3. **The check set is the main cost.** With 10 to 30 toxic compounds in the check set, the miss rate across repeated draws ranged from 0% to about 22% (5th to 95th percentile). Getting 50 toxic compounds into a random check set needs about 1,350 lab tests at 4% toxic and about 320 at 16% toxic.
4. **Stopping the first time the model "looks good" weakens the promise.** On SR-p53 with 30 toxic check compounds, the mean miss rose from 6.4% to 8.3%, and the share of runs missing more than 15% rose from 2% to 19%. The cautious alternative (3 pre-chosen looks, promise split 3 ways) was so conservative that it flagged 99% of compounds.
5. **New chemistry:** the similarity alarm fired on 100% of new-chemistry batches (false alarms 1% to 11%, nominal 5%). The miss rate did not rise (4% to 9%); the workload did (81% to 94% flagged, against 49% to 86% for familiar chemistry). Restarting, by recalibrating on 590 to 720 new lab tests, reduced the workload by 3 to 11 points.

## Limits (please read)
- Retrospective, on public data with a random forest and four endpoints. Not AstraZeneca's molecules.
- Tox21 labels stand in for "the lab"; assay noise and cost differences are not modelled.
- The simulated miss rate (about 6% to 7%) sits below the theoretical value (about 10%). Tied random forest scores may make the method conservative, but I have not tested that explanation.
- The plateau stopping rule in results/summary.md (section A2) is crude: it sometimes stops at a worse workload than running to the end. A better stopping rule is open work.
- The corrected stopping rule in experiment C is a conservative heuristic, not a proven guarantee.
- Novelty: I did not find a prior study combining all of these pieces, but my search was limited. Closest prior work: Corrigan et al. 2020 (conformal batch active learning), Bailey et al. (Sanofi batch active learning), Zhang, Norinder and Svensson 2021 (conformal toxicity), and a 2026 arXiv preprint on minority-class undercoverage.

## Run it
```
pip install -r requirements.txt
./get_data.sh
python run_experiments.py --exp A,B,C,D      # about 12 minutes on 2 cores; add --quick for a 1 minute test
python make_figures.py                       # writes figures/ and results/summary.md
```
Files: `wsal.py` (library), `run_experiments.py`, `make_figures.py`, `results/` (CSVs and summary), `figures/`.
