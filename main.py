import matplotlib.pyplot as plt
import numpy as np
from math import factorial
plt.style.use('_mpl-gallery')


def f(x):
    return 3.7 * np.cos((3 * x) / 2) * np.exp((-4 * x) / 3) + 2.4 * np.sin((9 * x) / 2) * np.exp((2 * x) / 3) + 4


def F(x):
    return f(x) / ((2.3 - x) ** (0.6))


def nodeEquidistant(leftBorder, rightBorder, nodeNumber):
    nodeList = []
    for i in range(0, nodeNumber):
        nodeList.append(leftBorder + i * (rightBorder - leftBorder) / (nodeNumber-1))
    return nodeList


def nodeOptimal(leftBorder, rightBorder, nodeNumber):
    nodeList = []
    for i in range(nodeNumber):
        nodeList.append(0.5 * ((-rightBorder + leftBorder) * np.cos(
            ((2 * i + 1) * np.pi) / (2 * nodeNumber)) + rightBorder + leftBorder))
    return nodeList


def integralWithSumDarboux(function, leftBorder, rightBorder, nodeDistribution, nodeNumber):
    result = 0
    nodeList = nodeDistribution(leftBorder, rightBorder, nodeNumber)
    for i in range(nodeNumber+1):
        result += function(nodeList[i]+(nodeList[i+1]-nodeList[i])/2)*(nodeList[i+1]-nodeList[i])
    return result


def printDarbouxResult(nodesNumber):
    integral = []
    for i in range(nodesNumber):
        integral.append(integralWithSumDarboux(F, 1.8, 2.3, nodeEquidistant, i))
    x = np.linspace(1.8, 2.2, nodesNumber)
    fig, ax = plt.subplots()
    ax.plot(x, integral, linewidth=1)
    plt.show()
    return integral[-1]


def polyCalc(n):
    polyWithCoef = []
    for i in range(n+1):
        #First element is coef, second is power
        coef = ((-1)**(i % 2))*(1/(i+2/5))*(2.3 ** (n - i)) * factorial(n)/(factorial(n-i)*factorial(i))
        polyWithCoef.append([coef, i+2/5])
    return polyWithCoef


def polyMean(poly, leftBorder, rightBorder):
    result = 0
    for monom in poly:
        result += -(monom[0]*((2.3-rightBorder) ** monom[1] - (2.3-leftBorder)**monom[1]))
    return result


def momentsCalculation(leftBorder, rightBorder, nodeNumber):
    moments = []
    for i in range(nodeNumber):
        poly = polyCalc(i)
        moments.append(polyMean(poly, leftBorder, rightBorder))
    return moments


def interpolationQuadratureFormula(nodeNumber, nodeDistribution):
    #Замена t=a-x следовательно в функцию подставлять сдвинутые значения
    result = 0
    moments = momentsCalculation(1.8, 2.3, nodeNumber)
    nodeList = nodeDistribution(1.8, 2.3, nodeNumber)
    X = np.transpose(np.vander(nodeList, increasing=True))
    A = np.linalg.solve(X, moments)
    for i in range(nodeNumber):
        result += A[i]*f(nodeList[i])
    return result


print(interpolationQuadratureFormula(20, nodeEquidistant))
#print(interpolationQuadratureFormula(3, nodeEquidistant))
#print(momentsCalculation(0,0.5, 1))git
#print(integralWithSumDarboux(F, 1.8, 2.3, nodeEquidistant, 500000))
