import matplotlib.pyplot as plt
import numpy as np
from math import factorial
import matplotlib as mpl

mpl.use('Qt5Agg')

integralCorrectMean = 1.185141974956241824914878594317090726677


def f(x):
    return 3.7 * np.cos((3 * x) / 2) * np.exp((-4 * x) / 3) + 2.4 * np.sin((9 * x) / 2) * np.exp((2 * x) / 3) + 4


def F(x):
    return f(x) / ((2.3 - x) ** 0.6)


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
    for i in range(nodeNumber-1):
        result += function(nodeList[i]+(nodeList[i+1]-nodeList[i])/2)*(nodeList[i+1]-nodeList[i])
    return result


def printDarbouxResult(nodesNumber):
    integral = []
    for i in range(3, nodesNumber):
        integral.append(integralWithSumDarboux(F, 1.8, 2.3, nodeEquidistant, i))
    x = np.linspace(3, nodesNumber, nodesNumber-3)
    plt.plot(x, integral, linewidth=1)
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


def interpolationQuadratureFormula(leftBorder, rightBorder, nodeNumber, nodeDistribution):
    result = 0
    moments = momentsCalculation(leftBorder, rightBorder, nodeNumber)
    nodeList = nodeDistribution(leftBorder, rightBorder, nodeNumber)
    X = np.transpose(np.vander(nodeList, increasing=True))
    A = np.linalg.solve(X, moments)
    for i in range(nodeNumber):
        result += A[i]*f(nodeList[i])
    return [result, A]


def quadratureFormulaGauss(leftBorder, rightBorder, nodeNumber):
    moments = momentsCalculation(leftBorder, rightBorder, 2*nodeNumber)
    right_part = [-moments[i] for i in range(nodeNumber, 2*nodeNumber)]
    left_part = []
    for i in range(nodeNumber):
        curr = []
        for j in range(nodeNumber):
            curr.append(moments[j+i])
        left_part.append(curr)
    a = np.linalg.solve(left_part, right_part)
    poly = [1]
    for i in range(nodeNumber):
        poly.append(a[nodeNumber-i-1])
    roots = np.roots(poly)
    for root in roots:
        if root < leftBorder or root > rightBorder:
            return None
    X = np.transpose(np.vander(roots, increasing=True))
    A = np.linalg.solve(X, moments[:nodeNumber])
    for a in A:
        if a < 0:
            return None
    result = 0
    for i in range(nodeNumber):
        result += A[i]*f(roots[i])
    return result, A


def printQF(nodesNumber):
    integral = []
    A_sum = []
    min = 100
    index = 0
    for i in range(3, nodesNumber):
        elem, A = interpolationQuadratureFormula(1.8, 2.3, i, nodeEquidistant)
        elem = np.log10(abs(elem-integralCorrectMean))
        integral.append(elem)
        A_sum.append(sum(abs(A)))
        if elem < min:
            min = elem
            index = i
    x = np.linspace(3, nodesNumber, nodesNumber - 3)
    fig, axs = plt.subplots(2)
    axs[0].plot(x, integral, linewidth=1)
    axs[0].set_title('Error rate')
    axs[1].plot(x, A_sum)
    axs[1].set_title('A sum')
    plt.show()
    return min, index


def printQFGauss(nodesNumber):
    integral = []
    A_sum = []
    min = 100
    index = 0
    for i in range(3, nodesNumber):
        elem, A = quadratureFormulaGauss(1.8, 2.3, i)
        elem = np.log10(abs(elem-integralCorrectMean))
        integral.append(elem)
        A_sum.append(sum(abs(A)))
        if elem < min:
            min = elem
            index = i
    x = np.linspace(3, nodesNumber, nodesNumber - 3)
    fig, axs = plt.subplots(2)
    axs[0].plot(x, integral, linewidth=1)
    axs[0].set_title('Error rate')
    axs[1].plot(x, A_sum)
    axs[1].set_title('A sum')
    plt.show()
    return min, index


def compoundQuadratureFormulas(nodeNumber, pointsPerSection, flag_is_newton):
    result = 0
    nodeList = nodeEquidistant(1.8, 2.3, nodeNumber)
    if flag_is_newton:
        for i in range(nodeNumber-1):
            result += interpolationQuadratureFormula(nodeList[i], nodeList[i+1], pointsPerSection, nodeEquidistant)
    else:
        for i in range(nodeNumber-1):
            result += quadratureFormulaGauss(nodeList[i], nodeList[i+1], pointsPerSection)
    return result


def processAitken(sum_h1, sum_h2, sum_h3, l):
    return -1*(np.log((sum_h3-sum_h2) / (sum_h2-sum_h1))) / (np.log(l))


def ruleRunge(m, l, sum_h1, sum_h2):
    r_h1 = (sum_h2 - sum_h1) / (1 - l ** (-m))
    r_h2 = (sum_h2 - sum_h1) / (l ** (-m) - 1)
    return [r_h1,r_h2]


def methodRichardson(sum_h1, sum_h2, l, m):
    return sum_h2 + ((sum_h2 - sum_h1) / (l**m - 1))


def optimalIntegral(h0, l, flag_is_newton):
    if flag_is_newton:
        pointsPerSection = 3
    else:
        pointsPerSection = 4
    flag = True
    prev_h1 = compoundQuadratureFormulas((0.5/h0), pointsPerSection, flag_is_newton)
    prev_h2 = compoundQuadratureFormulas((0.5/(h0*l)), pointsPerSection, flag_is_newton)
    prev_h3 = compoundQuadratureFormulas((0.5/(h0*l*l)), pointsPerSection, flag_is_newton)
    while flag:



        if :
            flag = False


#def printCompoundQF(nodesNumber):
#    integral = []
#    for i in range(3, nodesNumber):
#        integral.append(compoundQuadratureFormulas(i, nodeEquidistant, 3))
#    x = np.linspace(3, nodesNumber, nodesNumber - 3)
#    plt.plot(x, integral, linewidth=1)
#    plt.show()
#    return integral[-1]


#print(quadratureFormulaGauss(1.8, 2.3, 7)-integralCorrectMean)
#printCompoundQF(100)
#print(printQF(25))
print(printQFGauss(25))
#printDarbouxResult(10000)
#print(compoundQuadratureFormulas(10000, 3,)-integralCorrectMean)
#print(np.log10(interpolationQuadratureFormula(1.8,2.3,21, nodeEquidistant)-integralCorrectMean))
#print(interpolationQuadratureFormula(3, nodeEquidistant))
#print(momentsCalculation(0,0.5, 1))
#print(integralWithSumDarboux(F, 1.8, 2.3, nodeEquidistant, 1000000)-integralCorrectMean)

#Для составной квадратурной формуле на базе Ньютона котесса для отрезка берем 3, а для Гаусса 4