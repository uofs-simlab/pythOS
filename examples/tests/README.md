# Basic Testing Example for pythOS

This directory provides some basic tests of pythOS functionality to ensure the code and any desired dependencies are setup correctly.  Provided here is information about the options to control selection of dependencies and output.

## Fractional Step Tests

The file is `fractional_step_test.py`, which can be run as
```
python3 fractional_step_test.py [-v|--verbose] [-p|--plot] [-s|--sundials]
```

The options are:

`-v`, `--verbose`: Display the solutions at the end of solution time for all solution methods

`-p`, `--plot`: Display the plot of the time series for all selected solution methods

`-s`, `--sundials`: Test using the Sundials solvers with the adaptive solver and the EPI solver

The expected verbose output is:
```
adaptive solution                        [1.60065186 1.4930342 ]
basic split solution                     [1.66673022 1.4492337 ] with error [-0.06607837  0.0438005 ]
using fully implicit solver solution     [1.70234527 1.44546737] with error [-0.10169342  0.04756683]
using adaptive solver solution           [1.64991425 1.49525643] with error [-0.04926239 -0.00222223]
using Cash-Karp solver solution          [1.64992909 1.49524665] with error [-0.04927723 -0.00221245]
using CVODE solver solution              [1.64853688 1.49609869] with error [-0.04788503 -0.00306449]
using EPI (kiops) solver solution        [1.57627009 1.53343838] with error [ 0.02438177 -0.04040418]
using EPI (RK45) solver solution         [1.57627009 1.53343838] with error [ 0.02438177 -0.04040418]
using EPI (Cash-Karp) solver solution    [1.57627009 1.53343838] with error [ 0.02438177 -0.04040418]
using EPI (ARKODE) solver solution       [1.57627015 1.53343834] with error [ 0.0243817  -0.04040414]
basic (complex) split solution           [1.60114117-0.00137159j 1.49252805+0.00093171j] with error [-0.00048931+0.00137159j  0.00050615-0.00093171j]
using fully implicit solver solution     [1.62060617-0.00954216j 1.48170953+0.0053355j ] with error [-0.01995431+0.00954216j  0.01132467-0.0053355j ]
using Cash-Karp solver solution          [1.60057303+2.40076365e-05j 1.49304476-9.63393917e-07j] with error [ 7.88228383e-05-2.40076365e-05j -1.05562402e-05+9.63393917e-07j]
using EPI (kiops) solver solution        [1.60081695-0.00022341j 1.49291944+0.00013411j] with error [-0.00016509+0.00022341j  0.00011476-0.00013411j]
using EPI (RK45) solver solution         [1.60081695-0.00022341j 1.49291944+0.00013411j] with error [-0.00016509+0.00022341j  0.00011476-0.00013411j]
using EPI (Cash-Karp) solver solution    [1.60081695-0.00022341j 1.49291944+0.00013411j] with error [-0.00016509+0.00022341j  0.00011476-0.00013411j]
```

## Additive RK and Generalized Additive RK Tests

This script can be run as

```
python3 additive_test.py [-v|--verbose] [-p|--plot]
```

The options are:

`-v`, `--verbose`: Display the solutions at the end of solution time for all solution methods

`-p`, `--plot`: Display the plot of the time series for all selected solution methods

The expected verbose output is:

```
adaptive solution                        [1.60065186 1.4930342 ]
Additive RK solution                     [1.61722986 1.48315188] with error [-0.01657801  0.00988232]
Adaptive additive RK solution            [1.5969258  1.49532789] with error [ 0.00372605 -0.00229369]
Generalized additive RK solution         [1.74870444 1.40884735] with error [-0.14805258  0.08418685]
```

## Multirate and Multirate-Infinitesimal Tests

This script can be run as
```
python3 multirate_test.py [-v|--verbose] [-p|--plot] [-s|--sundials]
```

The options are:

`-v`, `--verbose`: Display the solutions at the end of solution time for all solution methods

`-p`, `--plot`: Display the plot of the time series for all selected solution methods

`-s`, `--sundials`: Test using the Sundials solver with the multirate infinitesimal methods

The expected verbose output is:
```
adaptive solution              [1.60065186 1.4930342 ]
mrgark_ex2_im2 solution        [1.60065186 1.4930342 ] with error [-0.1031802   0.06905785]
mri_kw3 solution               [1.60065186 1.4930342 ] with error [-0.04564358  0.03069739]
mri_kw3 (Cash-Karp) solution   [1.60065186 1.4930342 ] with error [-0.04564361  0.03069741]
mri_kw3 (CVode) solution       [1.60065186 1.4930342 ] with error [-0.04564361  0.03069741]
mri_imex3 solution             [1.60065186 1.4930342 ] with error [-0.0412165   0.04001363]
```
