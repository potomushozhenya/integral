import math
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
    if nodeNumber == 1:
        return [leftBorder+rightBorder/2]
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
    roots = np.sort(np.roots(poly))
    #for root in roots:
    #    if root < leftBorder or root > rightBorder:
    #        return None
    X = np.transpose(np.vander(roots, increasing=True))
    A = np.linalg.solve(X, moments[:nodeNumber])
    #for a in A:
    #    if a < 0:
    #        return None
    result = 0
    for i in range(nodeNumber):
        result += A[i]*f(roots[i])
    return [result, A]


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
            result += interpolationQuadratureFormula(nodeList[i], nodeList[i+1], pointsPerSection, nodeEquidistant)[0]
    else:
        for i in range(nodeNumber-1):
            result += quadratureFormulaGauss(nodeList[i], nodeList[i+1], pointsPerSection)[0]
    return result


def processAitken(sum_h1, sum_h2, sum_h3, l):
    return -1*(np.log((sum_h3-sum_h2) / (sum_h2-sum_h1))) / (np.log(l))


def ruleRunge(m, l, sum_h1, sum_h2, h_form, eps):
    return h_form / (0.95 * ( ((eps * (1 - l**(-m))) / abs(sum_h2 - sum_h1) )**(1/m)))


def optimalIntegral(h0, l, flag_is_newton, eps):
    #bigMesh булевая переменная отвечающая за сетки если true то h b h/l если false то h/l и h/l^2
    if flag_is_newton:
        pointsPerSection = 3
    else:
        pointsPerSection = 4
    h_form = math.ceil(0.5/h0)
    s = [compoundQuadratureFormulas(math.ceil(h_form * (l**i)), pointsPerSection, flag_is_newton) for i in range(3)]
    mPrev = 0
    while True:
        m = processAitken(s[0], s[1], s[2], l)
        print(h_form*l*l, m, integralCorrectMean - s[2])
        if abs(mPrev - m) < 0.15:
            break
        mPrev = m
        h_form = h_form * l
        s[0] = s[1]
        s[1] = s[2]
        s[2] = compoundQuadratureFormulas(h_form * l * l, pointsPerSection, flag_is_newton)
    h_opt = math.ceil(ruleRunge(m, l, s[0], s[1], h_form, eps))
    print(h_opt)
    s[2] = compoundQuadratureFormulas(h_opt, pointsPerSection, flag_is_newton)
    s[1] = compoundQuadratureFormulas(math.ceil(h_opt / l), pointsPerSection, flag_is_newton)
    s[0] = compoundQuadratureFormulas(math.ceil(h_opt / (l * l)), pointsPerSection, flag_is_newton)
    m = processAitken(s[0], s[1], s[2], l)
    h_opt = math.ceil(ruleRunge(m, l, s[0], s[1], h_opt / (l * l), eps))
    print(h_opt)
    s[2] = compoundQuadratureFormulas(h_opt, pointsPerSection, flag_is_newton)
    return s[2]
print(optimalIntegral(0.25, 2, 1, 1e-9) - integralCorrectMean)

#print(quadratureFormulaGauss(1.8, 2.3, 7)-integralCorrectMean)
#printCompoundQF(100)
#print(printQF(25))
#print(printQFGauss(25))
#printDarbouxResult(10000)
#print(compoundQuadratureFormulas(10000, 3,)-integralCorrectMean)
#print(np.log10(interpolationQuadratureFormula(1.8,2.3,21, nodeEquidistant)-integralCorrectMean))
#print(interpolationQuadratureFormula(3, nodeEquidistant))
#print(momentsCalculation(0,0.5, 1))
#print(integralWithSumDarboux(F, 1.8, 2.3, nodeEquidistant, 1000000)-integralCorrectMean)

#Для составной квадратурной формуле на базе Ньютона котесса для отрезка берем 3, а для Гаусса 4