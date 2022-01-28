import EigensetReader

e = EigensetReader.Eigenset()
e.read("results.eigenset")

for x in range(0, e.numberEigenvectors):
    print(str(e.eigenpairs[x].J) + " " + str(e.eigenpairs[x].Bx) + " " + str(e.eigenpairs[x].Bz) + " " + str(e.eigenpairs[x].Eigenvalue) )
