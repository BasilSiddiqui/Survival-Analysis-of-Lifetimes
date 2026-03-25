# 🐱 Survival Analysis of Cat Lifetimes (Censored Data)

This project explores survival patterns in two cat breeds using **time-to-event modelling, censoring techniques, and hazard-based inference**.

The analysis combines:
- Kaplan–Meier survival estimation  
- Conditional probability reasoning  
- Cox proportional hazards modelling  
- Illness-death modelling with infection dynamics  

---

## 📁 Data Overview

Two synthetic but realistic datasets were generated:

### 🐾 American Shorthair
- Sample size: ~160–180  
- Lifetime distribution: Gamma-based with censoring  
- Right-censored at 15 years  

### 🐾 British Shorthair
- Sample size: ~140–160  
- Slightly different lifetime distribution  
- Same censoring mechanism  

### ⚠️ Additional Variable
- **FIP infection time** (if applicable)
- Introduces a **state change**: healthy → infected → death  

---

## ⚙️ Methodology

### 1. Survival Analysis (Kaplan–Meier)
- Non-parametric estimation of survival functions
- Handles right-censored data
- Used to compute survival probabilities at specific times

### 2. Conditional Probability & Bayes Rule
- Probabilities derived from survival curves
- Bayesian reasoning applied to infer breed likelihood given death

### 3. Cox Proportional Hazards Model
- Estimates relative hazard between breeds
- Uses partial likelihood (no baseline hazard assumption)

### 4. Illness-Death Model
- Models transition:
  - Healthy → Infected → Death  
- Separate hazard rates:
  - μ (healthy death)
  - σ (infection rate)
  - ν (death after infection)

---

## 📊 Key Findings

### 🟢 Survival Probabilities

- **American cats surviving to 10 years**:  
  → ~12.5%  
  → Indicates relatively low long-term survival

- **British cats dying by 12 years**:  
  → ~83.5%  
  → Suggests slightly higher mortality over time

---

### 🔁 Conditional Survival

- Probability of death between ages 3 and 8:

| Breed | Probability |
|------|-----------|
| American | ~64.5% |
| British  | ~59.2% |

👉 American cats show slightly higher mid-life mortality.

---

### 🧠 Bayesian Inference

- Probability a cat is **American given it died by 8.2 years**:

→ **~54.9%**

👉 Despite similar populations, earlier mortality slightly increases the likelihood of being American.

---

### 📉 Cox Model Results

- Estimated coefficient (β): **-0.206**

👉 Interpretation:
- British cats have a **lower hazard rate** than American cats  
- (since β < 0)

---

### 🧪 Hypothesis Testing

- Likelihood Ratio Test statistic: **2.39**
- p-value: **0.122**

👉 Conclusion:
- No statistically significant difference in survival between breeds  
- (at 2% significance level)

---

### ⚠️ Model Limitation Insight

The Cox model assumes **constant hazard ratios over time**, but:

- Infection (FIP) changes risk dynamically  
- Hazard is **time-dependent**

👉 Therefore:
- Cox model is **not fully appropriate** in this setting

---

### 🦠 Infection Dynamics (Illness-Death Model)

- Estimated post-infection death rate (ν): **~0.99**

👉 Interpretation:
- Once infected, cats face **very high mortality risk**

---

### ⏱️ Infection Risk Over Time

- Probability of infection between ages 7 and 10:

→ **~90.4%**

👉 Indicates rapid increase in infection risk with age under this model.

---

## 🧠 Key Insights

- Survival patterns differ slightly between breeds, but not significantly  
- Mid-life mortality is high for both groups  
- Infection dramatically increases death risk  
- Ignoring time-dependent effects (like infection) can lead to misleading conclusions  
- Multi-state models provide a more realistic framework than standard survival models  

---

## 📈 Outputs

The script produces:
- Survival probabilities at key time points  
- Conditional probabilities  
- Cox model estimates and hypothesis tests  
- Illness-death likelihood estimates  
- Final results exported as csv
```

## ⚙️ Tech Used

- R  
- `survival` package  
- Base R for likelihood estimation and optimisation  

---

## 📌 Notes

- All data is reproducible via fixed random seed  
- No external dependencies beyond required libraries  
- Focus is on applied survival modelling and interpretation  

---

## 🚀 Takeaway

This project shows how survival analysis evolves when moving from:
- Simple time-to-event models  
→ to  
- State-dependent, real-world risk modelling  

Highlighting the importance of choosing the **right model for the data-generating process**.
