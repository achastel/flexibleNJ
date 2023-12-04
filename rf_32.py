from ete3 import Tree

NJPY=Tree("(((12,(((24,(((23,9),0),11)),(26,(30,15))),6)),7),((10,((((20,3),19),((17,21),22)),((25,2),((31,(14,29)),1)))),(((5,4),(27,28)),((18,13),(16,8)))));")

NJBIOR=Tree("((((((((30,15),26),(24,(((23,9),0),11))),6),12),7),(((28,27),(5,4)),((18,13),(16,8)))),(((((31,(29,14)),1),(25,2)),((22,(21,17)),((20,3),19))),10));")

NJBIOU=Tree("(((((((30,15),26),(24,(((23,9),0),11))),6),12),7),(((28,27),(5,4)),((18,13),(16,8))),(((((31,(29,14)),1),(25,2)),((22,(21,17)),((20,3),19))),10));")

FNJ_05=Tree("(((((0,13),(11,23)),((9,24),(15,30))),(((2,31),(6,26)),((7,12),(10,25)))),((((1,4),(14,29)),((16,18),(17,21))),(((3,22),(19,20)),((5,8),(27,28)))));")

FNJ_03=Tree("(0,((((1,(((2,31),12),((((9,24),(11,23)),((14,29),27)),(((15,30),26),((16,18),(17,21)))))),((((((3,22),(19,20)),25),(((5,8),28),13)),10),7)),4),6));")

FNJ_01=Tree("(0,((((((1,((3,(((9,24),((((((((((((11,23),(15,30)),((17,21),22)),26),((27,31),29)),25),(18,28)),19),20),16),13),14)),12)),10)),(5,8)),6),7),4),2));")

FNJ_005=Tree("((0,(1,(((((((((((3,(((((((((((((14,(29,31)),27),28),25),24),23),22),21),20),19),18),17),16)),(9,((15,30),26))),13),12),11),10),5),8),7),6),4))),2);")

FNJ_0025=Tree("((0,(1,(((((((((((3,(((((((((((((14,(29,31)),27),28),25),24),23),22),21),20),19),18),17),16)),(9,((15,30),26))),13),12),11),10),5),8),7),6),4))),2);")

FNJ_001=Tree("((0,(1,(((((((((((3,(((((((((((((14,(29,31)),27),28),25),24),23),22),21),20),19),18),17),16)),(9,((15,30),26))),13),12),11),10),5),8),7),6),4))),2);")

FNJ_0005=Tree("((0,(1,(((((((((((3,(((((((((((((14,(29,31)),27),28),25),24),23),22),21),20),19),18),17),16)),(9,((15,30),26))),13),12),11),10),5),8),7),6),4))),2);")

FNJSO_05=Tree("(((((12,7),(10,25)),((26,6),(31,2))),(((30,15),(24,9)),((23,11),(13,0)))),((((27,28),(8,5)),((20,19),(22,3))),(((21,17),(18,16)),((29,14),(4,1)))));")
FNJSO_03=Tree("(0,(6,(4,((7,(10,((13,(28,(8,5))),(25,((20,19),(22,3)))))),((((((21,17),(18,16)),(26,(30,15))),((27,(29,14)),((23,11),(24,9)))),(12,(31,2))),1)))));")
FNJSO_01=Tree("(0,(2,(4,(7,(6,((5,8),((10,((12,((14,(13,(16,(20,(19,((18,28),(25,((29,(31,27)),(26,((22,(21,17)),((30,15),(23,11)))))))))))),(24,9))),3)),1)))))));")
FNJSO_005=Tree("((((4,(6,(7,(8,(5,(10,(11,(12,(13,(((26,(30,15)),9),((16,(17,(18,(19,(20,(21,(22,(23,(24,(25,(28,(27,((31,29),14))))))))))))),3))))))))))),1),0),2);")
FNJSO_0025=Tree("((((4,(6,(7,(8,(5,(10,(11,(12,(13,(((26,(30,15)),9),((16,(17,(18,(19,(20,(21,(22,(23,(24,(25,(28,(27,((31,29),14))))))))))))),3))))))))))),1),0),2);")
FNJSO_001=Tree("((((4,(6,(7,(8,(5,(10,(11,(12,(13,(((26,(30,15)),9),((16,(17,(18,(19,(20,(21,(22,(23,(24,(25,(28,(27,((31,29),14))))))))))))),3))))))))))),1),0),2);")
FNJSO_0005=Tree("((((4,(6,(7,(8,(5,(10,(11,(12,(13,(((26,(30,15)),9),((16,(17,(18,(19,(20,(21,(22,(23,(24,(25,(28,(27,((31,29),14))))))))))))),3))))))))))),1),0),2);")

PY_BIOR = NJPY.robinson_foulds(NJBIOR)
PY_BIOU = NJPY.robinson_foulds(NJBIOU,unrooted_trees=True)
BIOR_U = NJBIOR.robinson_foulds(NJBIOU,unrooted_trees=True)

BIOR_PY = NJBIOR.robinson_foulds(NJPY)
BIOU_PY = NJBIOU.robinson_foulds(NJPY,unrooted_trees=True)
BIOU_R = NJBIOU.robinson_foulds(NJBIOR,unrooted_trees=True)

print("PY_BIOR ",PY_BIOR[0],"/",PY_BIOR[1])
print("PY_BIOU ",PY_BIOU[0],"/",PY_BIOU[1])
print("BIOR_U  ",BIOR_U[0],"/",BIOR_U[1])
print("BIOR_PY ",BIOR_PY[0],"/",BIOR_PY[1])
print("BIOU_PY ",BIOU_PY[0],"/",BIOU_PY[1])
print("BIOU_R  ",BIOU_R[0],"/",BIOU_R[1])

BIOU_FNJSO_05 = NJBIOU.robinson_foulds(FNJSO_05,unrooted_trees=True)
BIOU_FNJSO_03 = NJBIOU.robinson_foulds(FNJSO_03,unrooted_trees=True)
BIOU_FNJSO_01 = NJBIOU.robinson_foulds(FNJSO_01,unrooted_trees=True)
BIOU_FNJSO_005 = NJBIOU.robinson_foulds(FNJSO_005,unrooted_trees=True)
BIOU_FNJSO_0025 = NJBIOU.robinson_foulds(FNJSO_0025,unrooted_trees=True)
BIOU_FNJSO_001 = NJBIOU.robinson_foulds(FNJSO_001,unrooted_trees=True)
BIOU_FNJSO_0005 = NJBIOU.robinson_foulds(FNJSO_0005,unrooted_trees=True)
print("BIOU_FNJSO_05 ",BIOU_FNJSO_05[0],"/",BIOU_FNJSO_05[1])
print("BIOU_FNJSO_03 ",BIOU_FNJSO_03[0],"/",BIOU_FNJSO_03[1])
print("BIOU_FNJSO_01 ",BIOU_FNJSO_01[0],"/",BIOU_FNJSO_01[1])
print("BIOU_FNJSO_005 ",BIOU_FNJSO_005[0],"/",BIOU_FNJSO_005[1])
print("BIOU_FNJSO_0025 ",BIOU_FNJSO_0025[0],"/",BIOU_FNJSO_0025[1])
print("BIOU_FNJSO_001 ",BIOU_FNJSO_001[0],"/",BIOU_FNJSO_001[1])
print("BIOU_FNJSO_0005 ",BIOU_FNJSO_0005[0],"/",BIOU_FNJSO_0005[1])



BIOR_FNJSO_05 = NJBIOR.robinson_foulds(FNJSO_05)
BIOR_FNJSO_03 = NJBIOR.robinson_foulds(FNJSO_03)
BIOR_FNJSO_01 = NJBIOR.robinson_foulds(FNJSO_01)
BIOR_FNJSO_005 = NJBIOR.robinson_foulds(FNJSO_005)
BIOR_FNJSO_0025 = NJBIOR.robinson_foulds(FNJSO_0025)
BIOR_FNJSO_001 = NJBIOR.robinson_foulds(FNJSO_001)
BIOR_FNJSO_0005 = NJBIOR.robinson_foulds(FNJSO_0005)
print("BIOR_FNJSO_05 ",BIOR_FNJSO_05[0],"/",BIOR_FNJSO_05[1])
print("BIOR_FNJSO_03 ",BIOR_FNJSO_03[0],"/",BIOR_FNJSO_03[1])
print("BIOR_FNJSO_01 ",BIOR_FNJSO_01[0],"/",BIOR_FNJSO_01[1])
print("BIOR_FNJSO_005 ",BIOR_FNJSO_005[0],"/",BIOR_FNJSO_005[1])
print("BIOR_FNJSO_0025 ",BIOR_FNJSO_0025[0],"/",BIOR_FNJSO_0025[1])
print("BIOR_FNJSO_001 ",BIOR_FNJSO_001[0],"/",BIOR_FNJSO_001[1])
print("BIOR_FNJSO_0005 ",BIOR_FNJSO_0005[0],"/",BIOR_FNJSO_0005[1])
