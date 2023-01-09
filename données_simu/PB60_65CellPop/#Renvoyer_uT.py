#Renvoyer u+u_e
for i in range(0,101):

    fichier_u = "/Users/Flavien/Desktop/Pb65/TC_Pb65_"+str(i)+".vtk"
    fichier_ue = "/Users/Flavien/Desktop/Pb65/E_Pb65_"+str(i)+".vtk"

    fichier_uT = open("/Users/Flavien/Desktop/PB60_65CellPop/uT_Pb65_"+str(i)+".vtk", "w")


    f_u = open(fichier_u, "r")
    f_ue = open(fichier_ue, "r")
    line = 1
    limite = 5 + 6068 + 1 + 5889 +1 + 5889 + 3 + 1
    while line<limite:
        content_line = f_u.readline()
        f_ue.readline()
        fichier_uT.write(content_line)
        line = line +1
    while line <limite+5889:
        content_u = f_u.readline()
        content_ue = f_ue.readline()
        content_u.replace('E','*10**')
        content_ue.replace('E','*10**')
        content_u_float = float(content_u.replace('\n',''))
        content_ue_float = float(content_ue.replace('\n',''))
        fichier_uT.write(str(content_u_float+content_ue_float)+'\n')
        line = line + 1

    f_u.close()
    f_ue.close()
    fichier_uT.close()
