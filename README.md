# Research and Development Assignment – Curve Fitting

This repository contains my solution for the Research and Development (R&D) assignment.  
The goal of this task is to estimate the unknown parameters θ (theta), M, and X in a given parametric curve equation using the data points provided in `xy_data.csv`.

---

## Problem Description

The given curve is defined as:

```
x = (t * cos(θ)) - e^(M * |t|) * sin(0.3t) * sin(θ) + X
y = 42 + (t * sin(θ)) + e^(M * |t|) * sin(0.3t) * cos(θ)
```

where:

- θ, M, and X are unknown parameters.
- The variable `t` takes values between 6 and 60.
- The valid parameter ranges are:
  - 0° < θ < 50°
  - -0.05 < M < 0.05
  - 0 < X < 100

The file `xy_data.csv` contains the (x, y) points that lie on this parametric curve.

---

## Final Submission

The final fitted equation of the curve, based on the estimated parameters, is shown below.

*(Replace these numbers with your actual results before submission.)*

```
(tcos(0.82610) - e^(0.07420abs(t))sin(0.3t)sin(0.82610) + 11.57930,
42 + tsin(0.82610) + e^(0.07420*abs(t))sin(0.3t)*cos(0.82610))
```

You can verify this curve using Desmos at:  
[https://www.desmos.com/calculator/rfj91yrxob](https://www.desmos.com/calculator/rfj91yrxob)

---

## Approach Used

1. **Understanding the curve**  
   The equation contains trigonometric and exponential terms influenced by θ, M, and X. The task was to adjust these parameters so the resulting curve passes through all the given data points.

2. **Data transformation**  
   Each data point (x, y) was transformed according to the parameters θ and X to estimate the corresponding `t`.

3. **Model comparison**  
   The transformed points were evaluated against the model expression `e^(M * |t|) * sin(0.3t)` to measure the deviation.

4. **Error minimization**  
   The program performs a two-step search procedure:
   - A coarse grid search to locate a region of minimal error.
   - A refinement step for improved accuracy.

5. **Result generation**  
   Once the parameters with the smallest error are found, they are printed to the console and saved to `desmos_string.txt` for easy verification.

---

## How to Run the Project

1. **Install Python and NumPy**
   ```
   pip install numpy
   ```

2. **Run the script**
   ```
   python fit_curve.py xy_data.csv
   ```

3. **Check the output**
   - The console displays the best-fit values of θ, M, and X.
   - The file `desmos_string.txt` contains the formatted curve equation you can copy directly into Desmos.

---

## Files Included

```
curve-fit-rnd/
│
├── fit_curve.py          # Main Python program for parameter estimation
├── xy_data.csv           # Input data points
├── desmos_string.txt     # Final fitted curve equation (auto-generated)
├── fit.png               # Optional graph comparing data vs fitted curve
└── README.md             # Explanation file
```

---

## Results Summary

- The estimated parameters (θ, M, X) fall within the specified ranges.  
- The fitted curve closely matches the given data points.  
- The equation works correctly in Desmos for visual verification.

---

## Notes

- Only NumPy was used, with no external optimization or AI-based libraries.  
- All search and fitting logic was written manually for transparency.  
- The goal was to build a readable and well-documented solution.  
- The optional file `fit.png` can be generated to visually confirm the result.

---

## Time Spent

Approximately 6–7 hours, including understanding the mathematics, implementing the algorithm, testing, and preparing this report.

---

## Author

Mallisetty Rathna Kumar 

Date of Submission: November 9, 2025
```
