import EigensetReader

e = EigensetReader.Eigenset()
e.read("results.eigenset")

for x in range(0, e.numberEigenvectors):
    print(e.eigenpairs[x].J + " " + e.eigenpairs[x].Bx + " " + e.eigenpairs[x].Bz + " " + e.eigenpairs[x].Eigenvalue + "\n" )
