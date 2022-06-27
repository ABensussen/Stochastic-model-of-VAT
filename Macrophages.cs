// See https://aka.ms/new-console-template for more information
for (int k = 0; k <= 10; k++)
{
    string fileName = "Macrophages";
    string filepath = @"C:\Users\anton\OneDrive\Escritorio\" + fileName + "_" + k + "_.txt";
    string filepath1 = @"C:\Users\anton\OneDrive\Escritorio\" + fileName + "_" + k + "_pure_data_Stochastic.txt";

    //*****************
    #region
    //Step 1: Declaring control variables
    List<string> intialConditions = new List<string>();
    List<string> fixedPoints = new List<string>();
    List<string> fixedPoints1 = new List<string>();
    List<string> fixedPoints2 = new List<string>();
    List<string> fixedPoints3 = new List<string>();
    List<string> fixedPoints4 = new List<string>();

    StreamWriter outcomeModel = new StreamWriter(filepath);
    StreamWriter onlyDataModel = new StreamWriter(filepath1);

    List<string>[] fixe = { fixedPoints, fixedPoints1, fixedPoints2, fixedPoints3, fixedPoints4 };

    int maxN = 5;
    int maxN1 = 10000;
    int maxNtotal = (int)Math.Pow(2, maxN);
    int iterationN = 30;
    double noiseLimit = 0.08;
    string[] CellTags = { "M0", "M1", "M2", "M1-TAM", "M2-TAM" };

    string a1 = "0  0  0  0  0  0  0  0  1  1  0  1  0  0  1";
    string a2 = "1  0  0  0  1  1  1  1  1  0  0  0  0  0  0";
    string a3 = "0  0  0  0  0  0  0  0  1  1  1  1  1  0  1";
    string a4 = "1  0  0  0  1  0  0  0  1  0  1  1  1  1  0";
    string a5 = "0  0  0  0  1  0  0  0  1  1  1  1  1  0  1";

    intialConditions.Add(a1.Replace("  ", ""));
    intialConditions.Add(a2.Replace("  ", ""));
    intialConditions.Add(a3.Replace("  ", ""));
    intialConditions.Add(a4.Replace("  ", ""));
    intialConditions.Add(a5.Replace("  ", ""));

    #endregion
    //*****************

    //Step 3: Setting input for cases of study
    // Cellular inputs
    int IFNGe = 01;
    int IL12e = 01;
    int IFNBe = 0;
    int IL4e = 0;
    int TGFBe = 0;
    int IL10e = 0;
    int IL6e = 0; 
    int LPS = 0;
    int TNFe = 0;
    int IL1Be = 0;
    int GMCSF = 0;
    int Ceramide = 01;
    int Insulin = 01;


    for (int ft = 0; ft < maxN1; ft++)
    {
        //Step 4: Searching for the fixed points
        for (int h = 0; h < maxN; h++)
        {
            //*****************
            #region
            //Subprocess 1: Picking an initial condition
            string sentence = intialConditions[h];
            char[] charArr = sentence.ToCharArray();
            int[] xS = Array.ConvertAll(charArr, c => (int)Char.GetNumericValue(c));
            //Gene network
            List<int> IFNG = new List<int>();
            List<int> IFNB = new List<int>();
            List<int> IL12 = new List<int>();
            List<int> TBET = new List<int>();
            List<int> INOS = new List<int>();
            List<int> TNF = new List<int>();
            List<int> TLR4 = new List<int>();
            List<int> IL1B = new List<int>();
            List<int> IL6 = new List<int>();
            List<int> IL10 = new List<int>();
            List<int> IL4 = new List<int>();
            List<int> IRF4 = new List<int>();
            List<int> ARG1 = new List<int>();
            List<int> GATA3 = new List<int>();
            List<int> TGFB = new List<int>();

            IFNG.Add(xS[0]);
            IFNB.Add(xS[1]);
            IL12.Add(xS[2]);
            TBET.Add(xS[3]);
            INOS.Add(xS[4]);
            TNF.Add(xS[5]);
            TLR4.Add(xS[6]);
            IL1B.Add(xS[7]);
            IL6.Add(xS[8]);
            IL10.Add(xS[9]);
            IL4.Add(xS[10]);
            IRF4.Add(xS[11]);
            ARG1.Add(xS[12]);
            GATA3.Add(xS[13]);
            TGFB.Add(xS[14]);


            #endregion
            //****************

            for (int u = 1; u < iterationN; u++)
            {
                //***********
                double p0 = Pvalue(u);
                double p1 = Pvalue(u);
                double p2 = Pvalue(u);
                double p3 = Pvalue(u);
                double p4 = Pvalue(u);
                double p5 = Pvalue(u);
                double p6 = Pvalue(u);
                double p7 = Pvalue(u);
                double p8 = Pvalue(u);
                double p9 = Pvalue(u);
                double p10 = Pvalue(u);
                double p11 = Pvalue(u);
                double p12 = Pvalue(u);
                double p13 = Pvalue(u);
                double p14 = Pvalue(u);

                //Logic rule for IFNG
                if (IFNGe == 1 &&
                    TGFB[u - 1] == 0)
                {
                    if (p0 <= noiseLimit) { IFNG.Add(0); } else { IFNG.Add(1); }
                }
                else
                {
                    if (p0 <= noiseLimit) { IFNG.Add(1); } else { IFNG.Add(0); }

                }


                //Logic rule for IFNB
                if ((IFNBe == 1 ||
                    IFNB[u - 1] == 1) && !(
                    TGFB[u - 1] == 1 ||
                    IFNG[u - 1] == 1))
                {
                    if (p1 <= noiseLimit) { IFNB.Add(0); } else { IFNB.Add(1); }

                }
                else
                {
                    if (p1 <= noiseLimit) { IFNB.Add(1); } else { IFNB.Add(0); }

                }


                //Logic rule for GATA3
                if (TBET[u - 1] == 1 &&
                    GATA3[u - 1] == 0)
                {
                    if (p2 <= noiseLimit) { IL12.Add(0); } else { IL12.Add(1); }
                }
                else
                {
                    if (p2 <= noiseLimit) { IL12.Add(1); } else { IL12.Add(0); }

                }


                //Logic rule for TBET
                if (((IL12[u - 1] == 1 ||
                    IL12e == 1) && (
                    IFNG[u - 1] == 1 ||
                    IFNB[u - 1] == 1 ||
                    IL6[u - 1] == 1)) && !(
                    GATA3[u - 1] == 1 ||
                    TGFB[u - 1] == 1))
                {
                    if (p3 <= noiseLimit) { TBET.Add(0); } else { TBET.Add(1); }
                }
                else
                {
                    if (p3 <= noiseLimit) { TBET.Add(1); } else { TBET.Add(0); }
                }


                //Logic rule for iNOS
                if ((IFNG[u - 1] == 1 ||
                    IFNB[u - 1] == 1) && (
                    IL1B[u - 1] == 1 ||
                    TNF[u - 1] == 1 ||
                    IL6[u - 1] == 1 ||
                    TLR4[u - 1] == 1) || GMCSF == 1)
                {
                    if (p4 <= noiseLimit) { INOS.Add(0); } else { INOS.Add(1); }
                }
                else
                {
                    if (p4 <= noiseLimit) { INOS.Add(1); } else { INOS.Add(0); }
                }


                //Logic rule for TNF
                if ((TNFe == 1 || Insulin == 1 || Ceramide == 1 ||
                    TNF[u - 1] == 1 ||
                    TLR4[u - 1] == 1 ||
                    IL1B[u - 1] == 1 ||
                    IL6[u - 1] == 1) &&
                    IL10[u - 1] == 0 &&
                    IL4[u - 1] == 0)
                {
                    if (p5 <= noiseLimit) { TNF.Add(0); } else { TNF.Add(1); }
                }
                else
                {
                    if (p5 <= noiseLimit) { TNF.Add(1); } else { TNF.Add(0); }
                }


                //Logic rule for TLR4
                if ((LPS == 1 || Ceramide == 1) &&
                        IL4[u - 1] == 0)
                {
                    if (p6 <= noiseLimit) { TLR4.Add(0); } else { TLR4.Add(1); }
                }
                else
                {
                    if (p6 <= noiseLimit) { TLR4.Add(1); } else { TLR4.Add(0); }
                }


                //Logic rule for IL1B
                if ((IL1Be == 1 ||
                    IL1B[u - 1] == 1 ||
                    IL6[u - 1] == 1 ||
                    TNF[u - 1] == 1 ||
                    TLR4[u - 1] == 1) &&
                    IL10[u - 1] == 0 &&
                    IL4[u - 1] == 0)
                {
                    if (p7 <= noiseLimit) { IL1B.Add(0); } else { IL1B.Add(1); }
                }
                else
                {
                    if (p7 <= noiseLimit) { IL1B.Add(1); } else { IL1B.Add(0); }
                }


                //Logic rule for IL6
                if (IL6e == 1 ||
                    IL6[u - 1] == 1 ||
                    IL1B[u - 1] == 1 ||
                    TNF[u - 1] == 1 ||
                    TLR4[u - 1] == 1 ||
                    TGFB[u - 1] == 1)
                {
                    if (p8 <= noiseLimit) { IL6.Add(0); } else { IL6.Add(1); }
                }
                else
                {
                    if (p8 <= noiseLimit) { IL6.Add(1); } else { IL6.Add(0); }
                }


                //Logic rule for IL10
                if ((IL10e == 1 ||
                    IL10[u - 1] == 1 ||
                    IL6[u - 1] == 1 ||
                    TGFB[u - 1] == 1 ||
                    TLR4[u - 1] == 1) && !(
                    (IFNG[u - 1] == 1 ||
                    IFNB[u - 1] == 1)))
                {
                    if (p9 <= noiseLimit) { IL10.Add(0); } else { IL10.Add(1); }
                }
                else
                {
                    if (p9 <= noiseLimit) { IL10.Add(1); } else { IL10.Add(0); }
                }

                //Logic rule for IL4
                if ((IL4e == 1 ||
                    GATA3[u - 1] == 1) &&
                    TBET[u - 1] == 0)
                {
                    if (p10 <= noiseLimit) { IL4.Add(0); } else { IL4.Add(1); }
                }
                else
                {
                    if (p10 <= noiseLimit) { IL4.Add(1); } else { IL4.Add(0); }
                }

                //Logic rule for IRF4
                if ((IL4[u - 1] == 1 ||
                    IL6[u - 1] == 1 ||
                    IRF4[u - 1] == 1 || GMCSF == 1) && !(
                    TLR4[u - 1] == 1 ||
                    TNF[u - 1] == 1 ||
                    IL1B[u - 1] == 1))
                {
                    if (p11 <= noiseLimit) { IRF4.Add(0); } else { IRF4.Add(1); }
                }
                else
                {
                    if (p11 <= noiseLimit) { IRF4.Add(1); } else { IRF4.Add(0); }
                }

                //Logic rule for ARG1
                if (IRF4[u - 1] == 1 &&
                    IL4[u - 1] == 1)
                {
                    if (p12 <= noiseLimit) { ARG1.Add(0); } else { ARG1.Add(1); }
                }
                else
                {
                    if (p12 <= noiseLimit) { ARG1.Add(1); } else { ARG1.Add(0); }
                }

                //Logic rule for GATA3
                if ((GATA3[u - 1] == 1 ||
                    IL4[u - 1] == 1) && !(
                    TBET[u - 1] == 1 ||
                    TGFB[u - 1] == 1))
                {
                    if (p13 <= noiseLimit) { GATA3.Add(0); } else { GATA3.Add(1); }
                }
                else
                {
                    if (p13 <= noiseLimit) { GATA3.Add(1); } else { GATA3.Add(0); }
                }

                //Logic rule for TGFB
                if ((TGFBe == 1 ||
                    TGFB[u - 1] == 1 ||
                    IL4[u - 1] == 1 ||
                    IL6[u - 1] == 1 ||
                    IL10[u - 1] == 1) && !(
                    IFNG[u - 1] == 1 ||
                    IFNB[u - 1] == 1 ||
                    TLR4[u - 1] == 1 ||
                    TNF[u - 1] == 1))
                {
                    if (p14 <= noiseLimit) { TGFB.Add(0); } else { TGFB.Add(1); }
                }
                else
                {
                    if (p14 <= noiseLimit) { TGFB.Add(1); } else { TGFB.Add(0); }
                }
            }


            //Subprocess 3: Saving data
            var node1 = IFNG.Select(x => Convert.ToString(x)).ToList();
            var node2 = IFNB.Select(x => Convert.ToString(x)).ToList();
            var node3 = IL12.Select(x => Convert.ToString(x)).ToList();
            var node4 = TBET.Select(x => Convert.ToString(x)).ToList();
            var node5 = INOS.Select(x => Convert.ToString(x)).ToList();
            var node6 = TNF.Select(x => Convert.ToString(x)).ToList();
            var node7 = TLR4.Select(x => Convert.ToString(x)).ToList();
            var node8 = IL1B.Select(x => Convert.ToString(x)).ToList();
            var node9 = IL6.Select(x => Convert.ToString(x)).ToList();
            var node10 = IL10.Select(x => Convert.ToString(x)).ToList();
            var node11 = IL4.Select(x => Convert.ToString(x)).ToList();
            var node12 = IRF4.Select(x => Convert.ToString(x)).ToList();
            var node13 = ARG1.Select(x => Convert.ToString(x)).ToList();
            var node14 = GATA3.Select(x => Convert.ToString(x)).ToList();
            var node15 = TGFB.Select(x => Convert.ToString(x)).ToList();


            List<string> Nodes = new List<string>();

            for (int w = 0; w < node1.Count; w++)
            {
                Nodes.Add(node1[w] + "  " +
                                    node2[w] + "  " +
                                    node3[w] + "  " +
                                    node4[w] + "  " +
                                    node5[w] + "  " +
                                    node6[w] + "  " +
                                    node7[w] + "  " +
                                    node8[w] + "  " +
                                    node9[w] + "  " +
                                    node10[w] + "  " +
                                    node11[w] + "  " +
                                    node12[w] + "  " +
                                    node13[w] + "  " +
                                    node14[w] + "  " +
                                    node15[w]);

            }

            // Subprocess 4: Identifying an attractor (fixed points only)
            for (int z = 0; z < Nodes.Count - 1; z++)
            {
                if (Nodes[z + 1] == Nodes[z])
                {
                    fixe[h].Add(Nodes[z]);
                    break;
                }
            }

        }

    }

    for (int hn = 0; hn < maxN; hn++)
    {
        //Local
        int M0 = 0;
        int M1 = 0;
        int M2 = 0;
        int TAM_M1 = 0;
        int TAM_M2 = 0;
        int totalEnd = 0;
        //
        if (fixe[hn].Count == 0) { fixe[hn].Add("null"); }

        //Step 4: Calculating the size of the basin of attraction of each fixed point
        var q = from x in fixe[hn]
                group x by x into g
                let count = g.Count()
                orderby count descending
                select new { Value = g.Key, Count = count };

        //Step 5: Presentation of results
        outcomeModel.Write("      Nodes ID: {0  1  2  3  4  5  6  7  8  9  10 11 12 13 14} Initiated in  : {" + intialConditions[hn] + "}: " + CellTags[hn] + "\n");

        foreach (var x in q)
        {
            Console.WriteLine("Fixed point in: {" + x.Value + "}  Size of its basin of attraction: {" + x.Count + "}");
            outcomeModel.Write("Fixed point in: {" + x.Value + "}  Size of its basin of attraction: {" + x.Count + "}\n");
            //Labelling rules
            string ste = x.Value.Replace("  ", "");
            char[] Repase = ste.ToCharArray();
            Console.WriteLine(Repase.Length);

            if (Repase[4] == '0' &&
                Repase[12] == '0')
            {
                M0 = M0 + x.Count;
            }

            if (Repase[4] == '1' &&
                Repase[12] == '0')
            {
                M1 = M1 + x.Count;
            }

            if (Repase[4] == '0' &&
                Repase[12] == '1')
            {
                M2 = M2 + x.Count;
            }

            if (Repase[4] == '1' &&
                Repase[12] == '1' && (Repase[0] == '1' || Repase[2] == '1'))
            {
                TAM_M1 = TAM_M1 + x.Count;
            }

            if (Repase[4] == '1' &&
                Repase[12] == '1' && !(Repase[0] == '1' || Repase[2] == '1'))
            {
                TAM_M2 = TAM_M2 + x.Count;
            }

        }
        outcomeModel.Write("\n");


        //write

        outcomeModel.Write("\n");
        outcomeModel.Write("\n");
        outcomeModel.WriteLine("Total count: " + totalEnd);
        double resultM0 = 0;
        double resultM1 = 0;
        double resultM2 = 0;
        double resultTAMM1 = 0;
        double resultTAMM2 = 0;
        totalEnd = 0;
        totalEnd = M1 + M2 + M0 + TAM_M1 + TAM_M2;


        resultM0 = (double)M0 / (double)totalEnd;
        resultM1 = (double)M1 / (double)totalEnd;
        resultM2 = (double)M2 / (double)totalEnd;
        resultTAMM1 = (double)TAM_M1 / (double)totalEnd;
        resultTAMM2 = (double)TAM_M2 / (double)totalEnd;

        outcomeModel.Write("\n");

        outcomeModel.WriteLine("M0: " + M0);
        outcomeModel.WriteLine("M1: " + M1);
        outcomeModel.WriteLine("M2: " + M2);
        outcomeModel.WriteLine("TAM M1: " + TAM_M1);
        outcomeModel.WriteLine("TAM M2: " + TAM_M2);
        outcomeModel.Write("\n");

        outcomeModel.WriteLine(M0);
        outcomeModel.WriteLine(M1);
        outcomeModel.WriteLine(M2);
        outcomeModel.WriteLine(TAM_M1);
        outcomeModel.WriteLine(TAM_M2);


        onlyDataModel.WriteLine(resultM0);
        onlyDataModel.WriteLine(resultM1);
        onlyDataModel.WriteLine(resultM2);
        onlyDataModel.WriteLine(resultTAMM1);
        onlyDataModel.WriteLine(resultTAMM2);

        onlyDataModel.Write(" \n");
        onlyDataModel.Write(" \n");

    }

    outcomeModel.Write("\n");
    onlyDataModel.Write(" \n");
    onlyDataModel.Write(" \n");
    outcomeModel.Write("Total amount of initial states: " + maxNtotal.ToString());
    outcomeModel.Write("\n");
    outcomeModel.Write("Scalling factor: " + maxN1.ToString());
    outcomeModel.Write("\n");
    outcomeModel.Write("Legend:\n");
    outcomeModel.Write("\n");
    outcomeModel.Write(" IFNG: {0}\n IFNB: {1}\n IL12: {2}\n TBET: {3}\n INOS: {4}\n");
    outcomeModel.Write(" TNF: {5}\n TLR4: {6}\n IL1B: {7}\n IL6: {8}\n IL10: {9}\n");
    outcomeModel.Write(" IL4: {10}\n IRF4: {11}\n ARG1: {12}\n GATA3: {13}\n TGFB: {14}\n");

    outcomeModel.Write("\n");
    outcomeModel.Write("\n");
    outcomeModel.Close();
    onlyDataModel.Close();
    Console.WriteLine(" ");

}
Console.WriteLine("The report is ready");

static double Pvalue(int u)
{
    double p;
    Random rnd = new Random();

    if (u < 16)
    {
        p = rnd.NextDouble();
    }
    else
    {
        p = 1;
    }
    return p;
}

