#def printCompoundQF(nodesNumber):
#    integral = []
#    for i in range(3, nodesNumber):
#        integral.append(compoundQuadratureFormulas(i, nodeEquidistant, 3))
#    x = np.linspace(3, nodesNumber, nodesNumber - 3)
#    plt.plot(x, integral, linewidth=1)
#    plt.show()
#    return integral[-1]


#def methodRichardson(sum_h1, sum_h2, l, m):
#    return sum_h2 + ((sum_h2 - sum_h1) / (l**m - 1))