# https://github.com/shayfletcherz/NumericalAnalysisEx16.git

# EX16

import math
import sympy as sp
from sympy.utilities.lambdify import lambdify
from datetime import datetime

#Function for printing the original function and its derivative
def printDerived(function):
    x = sp.symbols('x')
    print("f(x) : ", function)
    derivedFunction = sp.diff(function, x)
    print("f'(x) : ", derivedFunction)
    secondDerived = sp.diff(derivedFunction)
    print("f''(x) : ", secondDerived)

#Derivative calculation
def calcDerived(function):
    x = sp.symbols('x')
    f_prime = function.diff(x)
    return f_prime

#The function returns the number of iterations to the crossing method
def calcError(startPoint, endPoint, epsilon):
    return -(math.log(epsilon / (endPoint - startPoint)) / math.log(2))

#The crossing method - will separate the ranges into smaller ranges and return the roots
def partition(polynom, startPoint, endPoint, method, choice, epsilon):
    day = str(datetime.today().day)
    hour = str(datetime.now().hour)
    minute = str(datetime.now().minute)
    x = sp.symbols('x')
    f = lambdify(x, polynom)
    polynomTag = calcDerived(polynom)
    fTag = lambdify(x, polynomTag)
    i = startPoint
    while i < endPoint: #Running on the ranges
        j = i + 0.1
        if choice == 1:  #Bisction
            error = calcError(startPoint, endPoint, epsilon)
            c, iteration = method(f, i, j, error, epsilon)
        elif choice == 2:  #Newton-Rapson
            c, iteration = method(polynom, i, j, epsilon)

        if c is not None:
            print("root: " + str(c) + "00000" + day + hour + minute + " , and " + str(iteration) + " iterations")

        #Actions on root
        if choice == 1:  #Bisction
            error = calcError(startPoint, endPoint, epsilon)
            c, iteration = method(fTag, i, j, error, epsilon)  # derived
        elif choice == 2:  #Newton-Rapson
            c, iteration = method(polynomTag, i, j, epsilon)

        if c is not None:
            if -epsilon < f(c) < epsilon:
                print("root: " + str(c) + "00000" + day + hour + minute + " , and " + str(iteration) + " iterations")
        i += 0.1

#Print function using Bisection
def bisectionMethodPrint(polynom, startPoint, endPoint, choice, epsilon):
    print("Bisection_Method\n")
    partition(polynom, startPoint, endPoint, runBisection, choice, epsilon)

#The function will display the amount of roots in the given range
def runBisection(function, startPoint, endPoint, error, epsilon):
    fXl = function(startPoint)
    fXr = function(endPoint)
    if (fXl * fXr) > 0:
        return None, None
    i = -1
    c = startPoint
    while endPoint - startPoint > epsilon and i < error:
        i += 1
        c = (startPoint + endPoint) / 2 # calculating middle point
        if (function(startPoint) * function(c)) > 0:
            startPoint = c
        else:
            endPoint = c
    if i + 1 > error: #If there are no roots
        print("Could not find root.")
        return None, None
    return c, i

#Printing function using Newton Rapson
def newtonRaphsonMethodPrint(polynom, startPoint, endPoint, choice, epsilon):
    print("Newton Raphson\n")
    partition(polynom, startPoint, endPoint, runNewtonRephson, choice, epsilon)

#The function will display the amount of roots in the given range
def runNewtonRephson(f, startA, endB, epsilon):
    x = sp.symbols('x')
    fTag = lambdify(x, calcDerived(f))
    f = lambdify(x, f)
    fXl = f(startA)
    fXr = f(endB)
    if (fXl * fXr) > 0:
        return None, None
    i = -1
    c = (startA + endB) / 2 #Midpoint calculation
    newC = abs(startA + endB)
    while i < 100:
        i += 1
        temp = newC
        newC = c - (f(c) / fTag(c))
        if abs(newC - c) < epsilon: #If there are no roots
            return newC, i
        c = temp
    print("Could not find root.")
    return None, None

#Function that returns area of the integral by simpson method
def simpson(function, startPoint, endPoint, parts):
    if parts % 2 == 1:  #check for even numbers of parts
        return None
    x = sp.symbols('x')
    func = lambdify(x, function)
    gap = abs(endPoint - startPoint) / parts  #h
    string = "Integral(" + str(startPoint) + ", " + str(endPoint) + ") = 1/3 * " + str(gap) + "[f(" + str(startPoint) + ")"
    appr = func(startPoint)
    for i in range(1, parts):
        if i % 2 == 0:  #
            string += " + 2 * f(" + str((i * gap) + startPoint) + ")"
            appr += 2 * func((i * gap) + startPoint)
        else:  # if not even
            string += " + 4 * f(" + str((i * gap) + startPoint) + ")"
            appr += 4 * func((i * gap) + startPoint)
        if i % 4 ==0:  #printing
            string += "\n"
    string += " * f(" + str(endPoint) + ")]\n"
    print(string)  #print
    appr += func(endPoint)
    appr *= 1 / 3 * gap
    return appr

#The function returns the area in the range by romberg method
def rombergMethod(function, startPoint, endPoint, limit, epsilon):
    results = [[0 for i in range(limit + 1)] for j in range(limit + 1)]  #creation of matrix
    for k in range(0, limit):
        res = trapezMethod(function, startPoint, endPoint, 2 ** k)  #calculate trapez method
        results[k+1][1] = res  #storing values
        print("R" + str(k+1) + "," + str(1) + " = " + str(res))  #print
    for j in range(2, limit + 1):
        for k in range(2, limit + 1):
            results[k][j] = results[k][j - 1] + ((1 / ((4 ** (j - 1)) - 1)) * (results[k][j - 1] - results[k - 1][j - 1]))
            print("R" + str(k) + "," + str(j) + " = " + str(results[k][j]))  #print
            if abs(results[k][j] - results[k - 1][j]) < epsilon:  #check if the difference is less then epsilon
                return results[k][j]

    print("\n Final result:\n")
    return results[j-1][k-1]

#function returns The area in the range by trapez method
def trapezMethod(function, startPoint, endPoint, segments):
    x = sp.symbols('x')
    function = lambdify(x, function)
    h = (endPoint - startPoint) / segments
    sum = 0
    while startPoint < endPoint: #run unti end point
        sum += 0.5 * ((startPoint + h) - startPoint) * (function(startPoint) + function(startPoint + h))
        startPoint += h
    return sum

def Main():
    day = str(datetime.today().day)
    hour = str(datetime.now().hour)
    minute = str(datetime.now().minute)
    x = sp.symbols('x')
    f = (x ** 2 * sp.exp(-x ** 2 + (5 * x) - 3)) * (3*x - 5)
    startPoint = 0
    endPoint = 3
    epsilon = 1.0e-12
    simpsonStart = 0.5
    simpsonEnd = 1

    printDerived(f)
    print("\n")
    bisectionMethodPrint(f, startPoint, endPoint, 1, epsilon)
    print("\n")
    newtonRaphsonMethodPrint(f, startPoint, endPoint, 2, epsilon)

    print("\n Simpson:\n")
    print(str(simpson(f, simpsonStart, simpsonEnd, 4)) + "00000" + day + hour + minute)
    print("\n Romberg:\n")
    print(str(rombergMethod(f, simpsonStart, simpsonEnd, 4, epsilon)) + "00000" + day + hour + minute)


Main()  #Run
