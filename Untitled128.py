#!/usr/bin/env python
# coding: utf-8

# # question 01
Write a Python function that takes in two arrays of data and calculates the F-value for a variance ratio
test. The function should return the F-value and the corresponding p-value for the test.
# In[1]:


import numpy as np
from scipy.stats import f

def variance_ratio_test(arr1, arr2):
    """
    Calculates the F-value and p-value for a variance ratio test between two arrays of data.
    
    Parameters:
    arr1 (array-like): The first array of data.
    arr2 (array-like): The second array of data.
    
    Returns:
    F (float): The F-value for the test.
    p (float): The p-value for the test.
    """
    # Calculate the variances of the two arrays
    var1 = np.var(arr1, ddof=1)
    var2 = np.var(arr2, ddof=1)
    
    # Ensure var1 is larger
    if var2 > var1:
        var1, var2 = var2, var1
        arr1, arr2 = arr2, arr1
    
    # Calculate the F-value and p-value
    F = var1 / var2
    df1 = len(arr1) - 1
    df2 = len(arr2) - 1
    p = f.sf(F, df1, df2)
    
    return F, p


# # question 02

# In[2]:


from scipy.stats import f

def critical_f(alpha, dfn, dfd):
    """
    Calculates the critical F-value for a two-tailed test given a significance level (alpha) and degrees of
    freedom for the numerator (dfn) and denominator (dfd) of an F-distribution.
    
    Parameters:
    alpha (float): The significance level.
    dfn (int): The degrees of freedom for the numerator.
    dfd (int): The degrees of freedom for the denominator.
    
    Returns:
    f_crit (float): The critical F-value.
    """
    f_crit = f.ppf(1 - alpha/2, dfn, dfd)
    
    return f_crit


# # question 03

# In[3]:


import numpy as np
from scipy.stats import f

# Set random seed for reproducibility
np.random.seed(123)

# Generate random samples from two normal distributions with known variances
mu1, mu2 = 0, 0
sigma1, sigma2 = 2, 2
n1, n2 = 30, 30

data1 = np.random.normal(mu1, sigma1, n1)
data2 = np.random.normal(mu2, sigma2, n2)

# Perform F-test to test for equal variances
F = np.var(data1, ddof=1) / np.var(data2, ddof=1)
dfn = n1 - 1
dfd = n2 - 1
p = f.sf(F, dfn, dfd) * 2

# Print results
print("F-value:", F)
print("Degrees of freedom:", dfn, ",", dfd)
print("p-value:", p)


# # question 04

# In[4]:


import numpy as np
from scipy.stats import f

# Set known variances and sample sizes
sigma1_sq = 10
sigma2_sq = 15
n1 = 12
n2 = 12

# Generate random samples from two normal distributions with known variances
mu1, mu2 = 0, 0
data1 = np.random.normal(mu1, np.sqrt(sigma1_sq), n1)
data2 = np.random.normal(mu2, np.sqrt(sigma2_sq), n2)

# Calculate sample variances and F-test statistic
s1_sq = np.var(data1, ddof=1)
s2_sq = np.var(data2, ddof=1)
F = s1_sq / s2_sq

# Calculate p-value and critical F-value
dfn = n1 - 1
dfd = n2 - 1
p_value = f.sf(F, dfn, dfd) * 2
f_crit = f.ppf(0.025, dfn, dfd)

# Compare F-test statistic to critical F-value and print results
if F > f_crit:
    print("The variances are significantly different (reject H0) with F-value of", F, "and p-value of", p_value)
else:
    print("The variances are not significantly different (fail to reject H0) with F-value of", F, "and p-value of", p_value)


# # question 05

# In[5]:


from scipy.stats import f

# Set claimed variance and sample size
sigma_sq = 0.005
n = 25

# Calculate sample variance and F-test statistic
s_sq = 0.006
F = s_sq / sigma_sq

# Calculate p-value and critical F-value
dfn = n - 1
dfd = float('inf')
p_value = f.sf(F, dfn, dfd)
f_crit = f.ppf(0.01, dfn, dfd)

# Compare F-test statistic to critical F-value and print results
if F > f_crit:
    print("The claim of variance being", sigma_sq, "is not justified (reject H0) with F-value of", F, "and p-value of", p_value)
else:
    print("The claim of variance being", sigma_sq, "is justified (fail to reject H0) with F-value of", F, "and p-value of", p_value)


# # question 06

# In[6]:


def f_distribution_mean_variance(dfn, dfd):
    """
    Calculates the mean and variance of an F-distribution given degrees of freedom for the numerator and denominator.
    Returns the mean and variance as a tuple.
    """
    if dfn <= 2:
        mean = None
    else:
        mean = dfd / (dfd - 2)
        
    if dfn <= 4:
        variance = None
    else:
        variance = (2 * dfd ** 2 * (dfn + dfd - 2)) / (dfn * (dfd - 2) ** 2 * (dfd - 4))
        
    return mean, variance


# # question 07
Q7. A random sample of 10 measurements is taken from a normal population with unknown variance. The
sample variance is found to be 25. Another random sample of 15 measurements is taken from another
normal population with unknown variance, and the sample variance is found to be 20. Conduct an F-test
at the 10% significance level to determine if the variances are significantly different.
To test if the variances of the two populations are significantly different, we can use an F-test. The null hypothesis is that the variances are equal, while the alternative hypothesis is that they are not equal.

The F-statistic can be calculated as the ratio of the two sample variances, with the larger variance in the numerator.

F = s1^2 / s2^2

Under the null hypothesis, this statistic follows an F-distribution with degrees of freedom (df1 = n1 - 1) and (df2 = n2 - 1), where n1 and n2 are the sample sizes of the two populations.

To conduct the F-test at the 10% significance level, we first need to find the critical F-value from the F-distribution with (df1 = 9) and (df2 = 14) degrees of freedom, and a significance level of 0.10.

Using a Python function to calculate the critical F-value, we get:
# In[7]:


from scipy.stats import f

def critical_f_value(alpha, dfn, dfd):
    """
    Returns the critical F-value from the F-distribution for a given significance level and degrees of freedom
    for the numerator and denominator.
    """
    return f.ppf(q=1-alpha/2, dfn=dfn, dfd=dfd)

alpha = 0.10
dfn1 = 9
dfd1 = 14

critical_f = critical_f_value(alpha, dfn1, dfd1)
print("Critical F-value:", critical_f)


# In[8]:


import numpy as np

s1_squared = 25
s2_squared = 20

F = s1_squared / s2_squared

print("F-value:", F)


# # question 08

# In[9]:


import numpy as np

restaurant_A = [24, 25, 28, 23, 22, 20, 27]
restaurant_B = [31, 33, 35, 30, 32, 36]

squared_differences_A = np.square(np.array(restaurant_A) - np.mean(restaurant_A))
squared_differences_B = np.square(np.array(restaurant_B) - np.mean(restaurant_B))

variance_A = np.sum(squared_differences_A) / (len(restaurant_A) - 1)
variance_B = np.sum(squared_differences_B) / (len(restaurant_B) - 1)

print("Sample variance of Restaurant A:", variance_A)
print("Sample variance of Restaurant B:", variance_B)


# In[10]:


F = variance_A / variance_B

print("F-value:", F)


# In[11]:


from scipy.stats import f

alpha = 0.05
dfn = len(restaurant_A) - 1
dfd = len(restaurant_B) - 1

critical_F = f.ppf(q=1-alpha/2, dfn=dfn, dfd=dfd)

print("Critical F-value:", critical_F)


# # question 09
Q9. The following data represent the test scores of two groups of students: Group A: 80, 85, 90, 92, 87, 83;
Group B: 75, 78, 82, 79, 81, 84. Conduct an F-test at the 1% significance level to determine if the variances
are significantly different.To conduct an F-test to determine if the variances of two groups of data are significantly different, we need to perform the following steps:

Step 1: Calculate the sample variances of the two groups.

Step 2: Calculate the F-statistic using the formula F = s1^2/s2^2, where s1^2 is the sample variance of Group A and s2^2 is the sample variance of Group B.

Step 3: Determine the critical F-value at the chosen level of significance and with degrees of freedom (df) equal to n1-1 and n2-1, where n1 and n2 are the sample sizes of Group A and Group B, respectively.

Step 4: Compare the calculated F-statistic to the critical F-value. If the calculated F-statistic is greater than the critical F-value, we reject the null hypothesis that the variances are equal; otherwise, we fail to reject the null hypothesis.

Let's perform these steps for the given data:

Step 1: Calculate the sample variances of the two groups.

Sample variance of Group A:
s1^2 = ((80-87.17)^2 + (85-87.17)^2 + (90-87.17)^2 + (92-87.17)^2 + (87-87.17)^2 + (83-87.17)^2) / (6-1) = 19.30

Sample variance of Group B:
s2^2 = ((75-80)^2 + (78-80)^2 + (82-80)^2 + (79-80)^2 + (81-80)^2 + (84-80)^2) / (6-1) = 8.00

Step 2: Calculate the F-statistic using the formula F = s1^2/s2^2.

F = 19.30/8.00 = 2.41

Step 3: Determine the critical F-value at the 1% significance level and with df = 5-1 = 4 and df = 6-1 = 5.

Using an F-table or calculator, we find that the critical F-value with df1 = 4 and df2 = 5 and alpha = 0.01 is 7.71.

Step 4: Compare the calculated F-statistic to the critical F-value.

Since the calculated F-statistic (2.41) is less than the critical F-value (7.71), we fail to reject the null hypothesis that the variances of the two groups are equal at the 1% significance level.

Therefore, we can conclude that there is no significant difference in the variances of the two groups.