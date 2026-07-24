# Colab to GitHub — walkthrough

One rule underneath all of it: **GitHub is where the work lives, Colab is a
rented machine.** Colab disconnects without warning and anything unpushed
is gone.

Do steps 1–3 once. Steps 4–9 are every session.

---

## 1. Put the repo on GitHub

On your own machine, in the folder containing `README.md`:

```bash
git init -b main
git add .
git commit -m "Reading notes, experiment 01 protocol, analysis notebook"
```

Create an **empty** repository at github.com/new — public, no README, no
.gitignore, no license, since you already have all three. Then:

```bash
git remote add origin https://github.com/<you>/protein-foundation-models.git
git push -u origin main
```

If you have the `gh` CLI, `gh repo create protein-foundation-models --public
--source=. --push` replaces the last two commands.

No local git? Use github.com/new, then "uploading an existing file", and drag
the whole folder in. Slower, works fine, loses nothing.

**Check:** the repo page shows `models/`, `experiments/`, `notebooks/`.

---

## 2. Notebook format — one strategy, not two

```bash
pip install jupytext
```

`.py` in jupytext percent format is canonical. `.ipynb` is a generated
artefact and is gitignored. That is the whole policy.

Earlier drafts also installed `nbstripout` and shipped a `.gitattributes`
filter for `*.ipynb`. That was dead configuration: a filter cannot run on
files git ignores. Both have been removed. One fewer local git filter is one
fewer "works on my machine" failure.

A notebook with outputs turns a one-word edit into a thousand-line diff, and
a diff nobody can read is a diff nobody can challenge. Keeping `.ipynb` out
of the repo achieves that without a filter.

---

## 3. Create a token and give it to Colab

Go to **github.com/settings/personal-access-tokens/new** (fine-grained, not
classic).

- Name: `colab-pfm`
- Expiration: 90 days
- Repository access: **Only select repositories** → this repo
- Permissions → Repository permissions → **Contents: Read and write**
- Nothing else. Metadata read-only is added automatically.

Generate, and copy it now — GitHub shows it once.

In Colab, click the **key icon** in the left sidebar → **Add new secret**:

- Name: `GH_TOKEN`
- Value: paste the token
- **Toggle "Notebook access" on.** This is the step people miss; without it
  `userdata.get` raises.

Never paste a token into a cell. Cell contents get saved and can get pushed.

---

## 4. Open a notebook and set the runtime

New Colab notebook. **Runtime → Change runtime type → CPU.** Experiment 01a
needs no GPU and burns quota for nothing if you leave it on T4.

---

## 5. First cell — authenticate and clone

```python
import pathlib, os
from google.colab import userdata

USER = "<your-github-username>"
REPO = "protein-foundation-models"

# Keeps the token out of clone URLs, cell output, and .git/config.
# It is still plaintext in /root/.git-credentials for the life of the
# runtime. On a single-user ephemeral VM the chmod is close to symbolic;
# the real protection is the narrow fine-grained token scope from step 3.
CREDS = pathlib.Path("/root/.git-credentials")
CREDS.write_text(f"https://x-access-token:{userdata.get('GH_TOKEN')}@github.com\n")
CREDS.chmod(0o600)
!git config --global credential.helper store
!git config --global user.name "<Your Name>"
!git config --global user.email "<your-email>"

if not pathlib.Path(REPO).exists():
    !git clone https://github.com/{USER}/{REPO}.git
%cd {REPO}
!git pull --rebase
!pip install -q -r requirements.txt
```

**Check:** `!ls` shows the repo contents and `!git status` says clean.

---

## 6. Run the experiment

```python
!python notebooks/01a_msa_depth_published_scores.py
```

Or work interactively:

```python
!jupytext --to notebook notebooks/01a_msa_depth_published_scores.py
```

then open the generated `.ipynb` from the file browser. It is gitignored, so
edit the `.py` and regenerate rather than the other way round.

**Check:** `results/esm2_ladder_by_depth.csv` and
`figures/esm2_scaling_by_depth.png` have new timestamps.

---

## 7. Provenance — already automatic

The script fetches from a pinned commit SHA and stamps that same variable, so
`results/provenance_01a_msa_depth.json` cannot drift from the input:

```python
PROTEINGYM_COMMIT = os.environ.get("PROTEINGYM_COMMIT", "144fe22b...")
assert PROTEINGYM_COMMIT != "main"
RAW = f"https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/{PROTEINGYM_COMMIT}"
```

To rerun against a newer ProteinGym, override rather than edit:

```python
%env PROTEINGYM_COMMIT=<new 40-char sha>
```

A provenance record has to describe what the code consumed, not what you
meant it to consume. Recording a SHA while fetching from a branch is the
failure this guards against.

---

## 8. Commit and push

Stage named files, not directories — a broad `git add results figures ...`
sweeps in unrelated work from other experiments:

```python
!git add \
    notebooks/01a_msa_depth_published_scores.py \
    results/esm2_ladder_by_depth.csv \
    results/esm2_15B_vs_650M_paired.csv \
    results/provenance_01a_msa_depth.json \
    figures/esm2_scaling_by_depth.png \
    experiments/01-msa-depth-confound/

!git diff --cached --stat     # review before committing
!git commit -m "01a: size x depth interaction, beta=-0.015 p=0.002 above 650M"
!git push
```

Messages carry the finding, not the action. "01a: interaction significant
above 650M" is useful in a year; "update notebook" is not.

**Check:** refresh the repo page. The commit is there and the figure renders.

---

## 9. Every session after this

```python
%cd protein-foundation-models
!git pull --rebase
```

first, push last. If you edited on your laptop between sessions and skip the
pull, you get a conflict at push time.

---

## Runtime limits

| Task | Runtime | Feasible |
|---|---|---|
| 01a, published scores | CPU | ~1 min |
| ESM-2 8M–650M scoring | T4 free | hours |
| ESM-2 3B | L4/A100 fp16 | Pro |
| ESM-2 15B | — | no |
| Evo 2 7B | — | no, custom CUDA kernels |

Free sessions die around 12 hours and often sooner. The rule is
**checkpoint, commit, push** — all three. A commit that never left the VM
dies with the VM, so a local commit protects against a failed cell and
nothing else.

```python
done = {p.stem for p in pathlib.Path("results/scores").glob("*.csv")}
for i, assay in enumerate(assays):
    if assay in done:
        continue
    ...
    df.to_csv(f"results/scores/{assay}.csv", index=False)   # survives a bad cell
    if i % 20 == 0:                                          # survives a dead VM
        !git add results/scores && git commit -q -m "01b checkpoint: through {assay}" || true
        !git push -q
```

Pushing per assay creates commit noise. A batch of 10–25 assays, or roughly
every 20–30 minutes, or one expensive model pass, are all reasonable units.

---

## What goes in

Commit: code, protocols, results CSVs, figures, `requirements.txt`.
Never: raw data, weights, `.ipynb`, tokens.

Results CSVs are small and belong in git. They let anyone check your figures
against your numbers without a GPU, which is the point.
