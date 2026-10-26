## Predictive Rule-Based Decision Tree for Agarwood Induction

This decision tree directly maps to the physical and chemical requirements a fungus needs to survive inside a living tree, allowing you to accurately classify unconfirmed fungal genomes based on their CAZyme profiles.

---

### Diagnostic Pathways

**Pathway A: The "Brute Force" Inducer**
These criteria ensure the fungus has the enzymatic arsenal to drill through solid, healthy wood and trigger massive chemical stress.

* **Total AA ≥ 130 & AA Rank Top 3:** *Aquilaria* wood is heavily fortified with lignin. The fungus must dedicate a massive percentage of its genome to Auxiliary Activities to blast through this armor.
* **Total GH ≥ 300:** Immense Glycoside Hydrolase firepower is needed to rapidly chop up cellulose once lignin is breached. This rapid destruction creates the chemical fragments that trigger the tree's resin defense.
* **Tier 1 Target ≥ 1:** Having at least one CAZyme Gene Cluster targeting lignin, xylan, or cellulose proves the fungus is evolutionarily specialized for structural wood rather than general soil decay.

**Pathway B: The "Wound Opportunist" Inducer**
These criteria identify fungi that cannot penetrate healthy bark but successfully parasitize pre-existing physical wounds.

* **Low Firepower Allowance:** Because the bark and lignin barrier is already physically broken, the fungus does not need high AA or GH counts to enter the tree.
* **Tier 2 Targets ≥ 3:** Inside a wound, the fungus feeds on living parenchyma cells and stored sap, requiring multiple clusters targeting pectin, arabinan, or starch. This parasitic feeding stresses the tree's metabolism.
* **Lack of Tier 1 Targets:** An absence of structural wood-decay clusters confirms the fungus is biologically restricted to the wound site and will not structurally rot the trunk.

**Pathway C: The Non-Inducer**
This pathway acts as a strict biological filter to remove false positives.

* **Non-Plant Substrates:** If a fungus lacks Tier 1 and Tier 2 targets but has high clusters for chitin or peptidoglycan, it primarily feeds on insects or microbes. Because it cannot digest plant tissue, it cannot trigger an agarwood response.

---

### Data-Driven Threshold Rationale

The numerical thresholds are derived directly from the biological dividing lines visible in the `confirmed_inducers_summary.tsv` data. This prevents the model from confusing a weak opportunist with an aggressive wood-rotter.

| Metric | Threshold | Biological Rationale & Data Reference |
| --- | --- | --- |
| **Total AA** | **≥ 130** | **The Armor Breaker Gap:** Places the cut-off exactly in the dead zone between the lowest true wood-decayer (*P. setosa* at 139) and the highest opportunist (*A. niger* at 118). Ensures a massive, dedicated investment in lignin destruction. |
| **Total GH** | **≥ 300** | **The Maceration Standard:** High-yield pathogens (*C. gloeosporioides*, *F. oxysporum*, *F. solani*, *B. dothidea*, *L. theobromae*) easily clear this hurdle, demonstrating the sheer enzymatic volume needed to rapidly dismantle healthy wood. |



