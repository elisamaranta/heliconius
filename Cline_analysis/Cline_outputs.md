# Cline analysis outputs

We want to model the residuals for corneal area (CA) and facet count (FC) by transect position. 

We use the HZAR  package in R and compare the fit for 5 different models using the R script Cline_analysis.R

model I: "none"A model with no exponential tails is desired
model II: "right"A model with just one exponential tail on the right is desired.
model III: "left"A model with just one exponential tail on the left is desired.
model IV: "mirror"A model with two exponential tails mirrored about the cline center is desired.
model V: "both"A model with two tails with independent parameters is desired.

## Facet count clines
>**log(Facet count) ~ log(Abdomen length) + sex residuals**

```R=
AICc
modelI   -1206.003
modelII  -1200.914
modelIII -1201.779
modelIV  -1201.892
modelV   -1197.488
```
Plotted together:
<img src="https://hackmd.io/_uploads/H1mN9JmDR.png" alt="image" />
![image](https://hackmd.io/_uploads/H1mN9JmDR.png)

>model I variation of parameters
```R=
> print(hzar.getLLCutParam(erato$analysis$model.selected,
+                          names(erato$analysis$model.selected$data.param)));
center2LLLow center2LLHigh width2LLLow width2LLHigh  muL2LLLow  muL2LLHigh   muR2LLLow muR2LLHigh   varL2LLLow varL2LLHigh
1  0.006018738      47.10231    15.49531     84.99815 -0.1420105 -0.02490493 0.003521054  0.0433219 4.503289e-08 0.002897112
   varR2LLLow varR2LLHigh   varH2LLLow varH2LLHigh
1 0.002381964 0.004820569 3.098752e-08 0.001932423

> print(hzar.get.ML.cline(erato$analysis$model.selected))
$param.free
       center    width         muL        muR         varL        varR         varH
1787 19.53007 76.54552 -0.06244081 0.01942549 0.0003618244 0.003505655 0.0004597017

$param.all
$param.all$center
[1] 19.53007

$param.all$width
[1] 76.54552

$param.all$muL
[1] -0.06244081

$param.all$muR
[1] 0.01942549

$param.all$varL
[1] 0.0003618244

$param.all$varR
[1] 0.003505655

$param.all$varH
[1] 0.0004597017


$clineFunc
function (x) 
-0.0624408131536537 + (0.0194254896565197 - -0.0624408131536537) * 
    (1/(1 + exp(-4 * (x - 19.5300743634062)/76.5455185507292)))
<bytecode: 0x000001de9f1de460>

$logLike
[1] 610.1405

$isValid
[1] TRUE

attr(,"class")
[1] "hzar.cline"
```

Best model:
![image](https://hackmd.io/_uploads/HkVB31QPA.png)

## Corneal area clines
>**log(corneal area) ~ log(Abdomen length) + sex residuals**
```R 
              AICc
modelI   -977.4293
modelII  -973.3504
modelIII -973.2569
modelIV  -973.2720
modelV   -969.3587
```
Plotted all together:
![image](https://hackmd.io/_uploads/Hk2WWWmwR.png)


Best model is model I. Here is the summary:
```R
> print(hzar.getLLCutParam(erato$analysis$model.selected,
+                          names(erato$analysis$model.selected$data.param)));
  center2LLLow center2LLHigh width2LLLow width2LLHigh  muL2LLLow muL2LLHigh  muR2LLLow muR2LLHigh  varL2LLLow varL2LLHigh
1     5.028478      47.94281    31.69893     84.99691 -0.1694051 -0.0467224 0.01344352 0.07391125 7.67994e-08 0.003107253
   varR2LLLow varR2LLHigh   varH2LLLow varH2LLHigh
1 0.002294954 0.008042034 5.714559e-07 0.004773607
> 
> ## Print the maximum likelihood cline for the selected model
> print(hzar.get.ML.cline(erato$analysis$model.selected))
$param.free
       center    width         muL        muR         varL        varR        varH
8777 29.94218 60.09856 -0.07927126 0.03445206 5.937826e-05 0.004936347 0.002564637

$param.all
$param.all$center
[1] 29.94218

$param.all$width
[1] 60.09856

$param.all$muL
[1] -0.07927126

$param.all$muR
[1] 0.03445206

$param.all$varL
[1] 5.937826e-05

$param.all$varR
[1] 0.004936347

$param.all$varH
[1] 0.002564637


$clineFunc
function (x) 
-0.0792712615097134 + (0.034452064517789 - -0.0792712615097134) * 
    (1/(1 + exp(-4 * (x - 29.9421788085021)/60.0985560114067)))

$logLike
[1] 495.8536

$isValid
[1] TRUE

attr(,"class")
[1] "hzar.cline"
```
...and plot:
![image](https://hackmd.io/_uploads/Bks6l-XPR.png)


Now, facet count (red) and corneal area (black) clines overlapping to visualise better:

![image](https://hackmd.io/_uploads/r1NTpQmD0.png)

Here is the WntA phenotype cline for comparison:

![image](https://hackmd.io/_uploads/rJ8Rh7XPA.png)


## Clines using tibia residuals
**Now we are look at tibia residuals instead of abdomen length residuals:**

log10(facet count) ~ log10(tibia length) + sex


```R 
              AICc
modelI   -1024.946
modelII  -1039.630
modelIII -1020.759
modelIV  -1021.550
modelV   -1035.428

```
![image](https://hackmd.io/_uploads/Hk3KRXrwC.png)

The best model is model II:
```R 
> deltaAICc
[1] 4.201981

> print(hzar.getLLCutParam(erato$analysis$model.selected,
+                          names(erato$analysis$model.selected$data.param)));
  center2LLLow center2LLHigh width2LLLow width2LLHigh  muL2LLLow  muL2LLHigh   muR2LLLow muR2LLHigh   varL2LLLow varL2LLHigh
1  0.005901026      72.17163   0.4125585     84.98728 -0.1692917 -0.01881676 0.003187303 0.09500472 1.574984e-05 0.006299549
    varR2LLLow varR2LLHigh   varH2LLLow varH2LLHigh deltaR2LLLow deltaR2LLHigh   tauR2LLLow tauR2LLHigh
1 0.0003360885 0.004080819 3.614329e-08 0.005335424   0.05777461      84.97035 6.097551e-05   0.9997941
> 
> ## Print the maximum likelihood cline for the selected model
> print(hzar.get.ML.cline(erato$analysis$model.selected))
$param.free
       center    width         muL        muR        varL        varR         varH   deltaR      tauR
6243 28.69351 76.33322 -0.06466693 0.03651217 0.002785873 0.002772965 4.015305e-05 10.84312 0.3411642

$param.all
$param.all$center
[1] 28.69351

$param.all$width
[1] 76.33322

$param.all$muL
[1] -0.06466693

$param.all$muR
[1] 0.03651217

$param.all$varL
[1] 0.002785873

$param.all$varR
[1] 0.002772965

$param.all$varH
[1] 4.015305e-05

$param.all$deltaR
[1] 10.84312

$param.all$tauR
[1] 0.3411642


$clineFunc
function (x) 
ifelse(x > 28.6935110330629 + 10.8431169239813, 0.0365121654607288 + 
    (-0.0646669293198544 - 0.0365121654607288) * (exp(4 * 10.8431169239813/76.3332187930192 * 
        (0.341164226079546 - 1))/(1 + exp(-(4 * 10.8431169239813/76.3332187930192))) * 
        exp(-4 * (x - 28.6935110330629)/76.3332187930192 * 0.341164226079546)), 
    -0.0646669293198544 + (0.0365121654607288 - -0.0646669293198544) * 
        (1/(1 + exp(-4 * (x - 28.6935110330629)/76.3332187930192))))
<bytecode: 0x0000014176df4458>

$logLike
[1] 529.0768

$isValid
[1] TRUE

attr(,"class")
[1] "hzar.cline"

```
Here is model II fit:
![image](https://hackmd.io/_uploads/SyVyRQSvA.png)

**Now repeat with corneal area:**

```R AICc
modelI   -876.1672
modelII  -866.4264
modelIII -871.6335
modelIV  -873.6386
modelV   -867.5695
```
![image](https://hackmd.io/_uploads/BkjaXBrPR.png)


Model I is the best fit cline:

```R
deltaAICc
[1] 2.528582

> print(hzar.getLLCutParam(erato$analysis$model.selected,
+                          names(erato$analysis$model.selected$data.param)));
  center2LLLow center2LLHigh width2LLLow width2LLHigh  muL2LLLow  muL2LLHigh  muR2LLLow muR2LLHigh  varL2LLLow varL2LLHigh   varR2LLLow
1     24.53379      84.94561    36.11355     84.99895 -0.1159188 -0.02792497 0.03044156  0.2340074 0.001616771 0.005701129 5.415907e-05
  varR2LLHigh   varH2LLLow varH2LLHigh
1  0.01059483 5.228082e-08  0.00368168
> 
> ## Print the maximum likelihood cline for the selected model
> print(hzar.get.ML.cline(erato$analysis$model.selected))
$param.free
       center   width         muL        muR        varL        varR         varH
1883 49.92212 84.1259 -0.06600927 0.08094324 0.003834766 0.005396273 0.0002101592

$param.all
$param.all$center
[1] 49.92212

$param.all$width
[1] 84.1259

$param.all$muL
[1] -0.06600927

$param.all$muR
[1] 0.08094324

$param.all$varL
[1] 0.003834766

$param.all$varR
[1] 0.005396273

$param.all$varH
[1] 0.0002101592


$clineFunc
function (x) 
-0.0660092740898278 + (0.0809432387548846 - -0.0660092740898278) * 
    (1/(1 + exp(-4 * (x - 49.9221150708329)/84.1258984407555)))

$logLike
[1] 445.2455

$isValid
[1] TRUE

attr(,"class")
[1] "hzar.cline"
```
And plotted here: 
![image](https://hackmd.io/_uploads/ryaKoBSwR.png)


Now together...

![image](https://hackmd.io/_uploads/rJr_iHBD0.png)


## Wing aspect ratio clines

I also tested the cline of wing aspect ratio to assess how it compared to eye size clines:

This plot shows wing aspect ratio along sites:

![image](https://hackmd.io/_uploads/SyWLdR1tA.png)

```R
              AICc
modelI   -859.3483
modelII  -853.6912
modelIII -840.7712
modelIV  -840.8987
modelV   -836.4745
```

All models overlapping:
![image](https://hackmd.io/_uploads/Hyup04eYA.png)

Best model fit is model I:

```R 
  center2LLLow center2LLHigh width2LLLow width2LLHigh muL2LLLow muL2LLHigh muR2LLLow muR2LLHigh   varL2LLLow varL2LLHigh  varR2LLLow varR2LLHigh
1 0.0009634791      8.017403    2.993656     36.56332 -78.16112   2.090051  2.136439   2.148718 0.0001246408    3.388493 0.003101569 0.004143033
   varH2LLLow varH2LLHigh
1 0.000104517   0.9217768

$param.free
      center    width      muL      muR        varL        varR        varH
8664 3.17972 10.78117 1.842877 2.144123 0.008262443 0.003465689 0.001356385

$param.all
$param.all$center
[1] 3.17972

$param.all$width
[1] 10.78117

$param.all$muL
[1] 1.842877

$param.all$muR
[1] 2.144123

$param.all$varL
[1] 0.008262443

$param.all$varR
[1] 0.003465689

$param.all$varH
[1] 0.001356385


$clineFunc
function (x) 
1.84287711051459 + (2.14412321965785 - 1.84287711051459) * (1/(1 + 
    exp(-4 * (x - 3.17971954003392)/10.7811667873813)))
<bytecode: 0x00000206a3444d48>

$logLike
[1] 436.8565

$isValid
[1] TRUE

attr(,"class")
[1] "hzar.cline"
```

![image](https://hackmd.io/_uploads/rynPkHxt0.png)
Not a good fit and no clear cline from the wing shape aspect ratio values. 


