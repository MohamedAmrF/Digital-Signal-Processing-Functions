# 1. Discrete Fourier Transform (DFT.cpp)
A program to compute DFT of a given signal.

## Input:
Enter the number of samples in your signal in the first line.  ( ***n*** )  
Enter the sample values of your signal in the second line.  ***X(0, 1, 2, .., n-1)***

### Example: 
```
4
1 1 -1 0 
```
## Output:
The output will contain the same number of samples after applying **DFT** to it.  
Each sample will be represented in complex numbers.
complex number will be written in the format of (real, imag).
### Example:
- 1 + 1i &rarr; (1, -1)
- -i     &rarr; (0, -1)

## Demo:
### input:
```
4
1 1 -1 0 
```
### output:
```
(1,0) (2,-1) (-1,0) (2,1) 
```

---------------------------

<br/>

# 2.Inverse Discrete fourier transform (IDFT.cpp)
A program to compute ***IDFT*** of a given signal.

## Input:
Enter the number of samples (***N***) in your signal in the first line.  
Enter the sample values of your signal in the following (***N***) lines in the following format: "real imag".

### Example: 
```
4
1 0 
2 -1
-1 0  
2 1 
```
## Output:
The output will contain the same number of samples after applying **IDFT** to it.  
Each sample will be represented in complex numbers.
complex number will be written in the format of (real, imag).
### Example:
- 1 + 1i &rarr; (1, -1)
- -i     &rarr; (0, -1)

## Demo:
### input:
```
4
1 0 
2 -1
-1 0  
2 1 
```
### output:
```
(1,0) (1,0) (-1,0) (0,0) 
```

<br> <br/>

# 3. Digital Convolution (convolution.cpp)

## Input:
Takes two inputs and convolves them.  
Two samples are represented in waveforms.c file.

## Output: 
Return a vector of samples representing the output signal.

Writes ***output signal*** to a ***.dat*** file.  
Writes ***input signal*** to a ***.dat*** file.  
Writes ***impulse response*** to a ***.dat*** file.  

- You can use **gnuplot** to plot the ***.dat*** files.  
  
<br> <br/>

# 4. Windowing Functions (windows.cpp)

- Available functions: ```triang```, ```hamming```, ``hanning``.  

## Input:
You Can construct a vector using the construct vector function.  
``` cpp
vector<double>t = construct_vector(initial_value, final_value, step);
vector<double>t = construct_vector(0, 10, 0.5);
```
this function simulates matlab's ```linspace``` function or this line in matlab: 
``` matlab
t = 0:0.5:10;
```

## Output: 
Prints the ***vector*** after using any window function.  
Creates files using ```write_to_file``` function as shown in the main as an example.    
Which can be plotted using ***gnuplot***.  