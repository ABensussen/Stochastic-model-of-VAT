// See https://aka.ms/new-console-template for more information
for (int k = 0; k <= 10; k++)
{
    string fileName = "Adipocytes";
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
    List<string> fixedPoints5 = new List<string>();
    List<string> fixedPoints6 = new List<string>();
    List<string> fixedPoints7 = new List<string>();

    StreamWriter outcomeModel = new StreamWriter(filepath);
    StreamWriter onlyDataModel = new StreamWriter(filepath1);

    List<string>[] fixe = { fixedPoints, fixedPoints1, fixedPoints2, fixedPoints3, fixedPoints4, fixedPoints5, fixedPoints6, fixedPoints7 };

    int maxN = 8;
    int maxN1 = 10000;
    int maxNtotal = (int)Math.Pow(2, maxN);
    int iterationN = 30;
    double noiseLimit = 0.08;
    string[] CellTags = {"TNF(-) CTGF(-) GLUT4(-)",
                                    "TNF(-) CTGF(-) GLUT4(+)",
                                    "TNF(-) CTGF(+) GLUT4(-)",
                                    "TNF(-) CTGF(+) GLUT4(+)",
                                    "TNF(+) CTGF(-) GLUT4(-)",
                                    "TNF(+) CTGF(-) GLUT4(+)",
                                    "TNF(+) CTGF(+) GLUT4(-)",
                                    "TNF(+) CTGF(+) GLUT4(+)" };

    string a1 = "0  0  0  0  1  1  0  1  0  0  0  1  0  0  0  0";
    string a2 = "0  0  0  0  1  1  0  1  0  1  0  1  0  0  0  1";
    string a3 = "0  0  0  0  1  1  0  1  1  0  0  1  0  0  0  0";
    string a4 = "0  0  0  0  1  1  0  1  1  1  0  1  0  0  0  1";
    string a5 = "1  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0";
    string a6 = "1  1  1  1  1  1  1  0  0  1  1  0  0  0  0  1";
    string a7 = "1  1  1  1  1  1  0  0  1  0  0  0  0  0  0  0";
    string a8 = "1  1  1  1  1  1  1  0  1  1  1  0  0  0  0  1";

    intialConditions.Add(a1.Replace("  ", ""));
    intialConditions.Add(a2.Replace("  ", ""));
    intialConditions.Add(a3.Replace("  ", ""));
    intialConditions.Add(a4.Replace("  ", ""));
    intialConditions.Add(a5.Replace("  ", ""));
    intialConditions.Add(a6.Replace("  ", ""));
    intialConditions.Add(a7.Replace("  ", ""));
    intialConditions.Add(a8.Replace("  ", ""));

    #endregion
    //*****************

    //Step 3: Setting input for cases of study
    // Cellular inputs
    int TNFe = 01;
    int TNFR2 = 01;
    int IL4e = 01;
    int IL6e = 01;
    int IL10e = 01;
    int IFNG = 01;
    int Insulin = 01;
    int CeramideExt = 01;
    int MechanicalPressure = 0;
    int CTGFR = 0;


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
            List<int> IL1b = new List<int>();
            List<int> IL6 = new List<int>();
            List<int> IL8 = new List<int>();
            List<int> TNF = new List<int>();
            List<int> IL10 = new List<int>();
            List<int> IL4 = new List<int>();
            List<int> AP1 = new List<int>();
            List<int> PPARg = new List<int>();
            List<int> CTGF = new List<int>();
            List<int> Akt = new List<int>();
            List<int> Resistin = new List<int>();
            List<int> Adiponectin = new List<int>();
            List<int> Ceramide = new List<int>();
            List<int> ROS = new List<int>();
            List<int> sFFA = new List<int>();
            List<int> GLUT4 = new List<int>();

            IL1b.Add(xS[0]);
            IL6.Add(xS[1]);
            IL8.Add(xS[2]);
            TNF.Add(xS[3]);
            IL10.Add(xS[4]);
            IL4.Add(xS[5]);
            AP1.Add(xS[6]);
            PPARg.Add(xS[7]);
            CTGF.Add(xS[8]);
            Akt.Add(xS[9]);
            Resistin.Add(xS[10]);
            Adiponectin.Add(xS[11]);
            Ceramide.Add(xS[12]);
            ROS.Add(xS[13]);
            sFFA.Add(xS[14]);
            GLUT4.Add(xS[15]);


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
                double p15 = Pvalue(u);

                //Logic rule for IL1B
                if ((Akt[u - 1] == 1 ||
                    Resistin[u - 1] == 1 ||
                    TNF[u - 1] == 1
                    ) && PPARg[u - 1] == 0)
                {
                    if (p0 <= noiseLimit) { IL1b.Add(0); } else { IL1b.Add(1); }
                }
                else
                {
                    if (p0 <= noiseLimit) { IL1b.Add(1); } else { IL1b.Add(0); }

                }


                //Logic rule for IL6
                if ((Akt[u - 1] == 1 ||
                    Resistin[u - 1] == 1 ||
                    TNF[u - 1] == 1 ||
                    IL6e == 1) && PPARg[u - 1] == 0)
                {
                    if (p1 <= noiseLimit) { IL6.Add(0); } else { IL6.Add(1); }

                }
                else
                {
                    if (p1 <= noiseLimit) { IL6.Add(1); } else { IL6.Add(0); }

                }


                //Logic rule for IL8
                if ((
                    TNF[u - 1] == 1 ||
                    Resistin[u - 1] == 1 ||
                    Akt[u - 1] == 1) && PPARg[u - 1] == 0)
                {
                    if (p2 <= noiseLimit) { IL8.Add(0); } else { IL8.Add(1); }
                }
                else
                {
                    if (p2 <= noiseLimit) { IL8.Add(1); } else { IL8.Add(0); }

                }


                //Logic rule for TNF
                if ((Akt[u - 1] == 1 ||
                    Resistin[u - 1] == 1 ||
                    (TNF[u - 1] == 1 && TNFR2 == 1) ||
                    (TNFe == 1 && TNFR2 == 1)) && PPARg[u - 1] == 0)
                {
                    if (p3 <= noiseLimit) { TNF.Add(0); } else { TNF.Add(1); }
                }
                else
                {
                    if (p3 <= noiseLimit) { TNF.Add(1); } else { TNF.Add(0); }
                }


                //Logic rule for IL10
                if ((Akt[u - 1] == 1 ||
                    IL10[u - 1] == 1 ||
                    IL4[u - 1] == 1 ||
                    IL10e == 1))
                {
                    if (p4 <= noiseLimit) { IL10.Add(0); } else { IL10.Add(1); }
                }
                else
                {
                    if (p4 <= noiseLimit) { IL10.Add(1); } else { IL10.Add(0); }
                }


                //Logic rule for IL4
                if ((IL4[u - 1] == 1 ||
                    IL10[u - 1] == 1 ||
                    IL4e == 1))
                {
                    if (p5 <= noiseLimit) { IL4.Add(0); } else { IL4.Add(1); }
                }
                else
                {
                    if (p5 <= noiseLimit) { IL4.Add(1); } else { IL4.Add(0); }
                }


                //Logic rule for AP1
                if ((Insulin == 1) && PPARg[u - 1] == 0)
                {
                    if (p6 <= noiseLimit) { AP1.Add(0); } else { AP1.Add(1); }
                }
                else
                {
                    if (p6 <= noiseLimit) { AP1.Add(1); } else { AP1.Add(0); }
                }


                //Logic rule for PPARg
                if ((IL4[u - 1] == 1 ||
                    IL10[u - 1] == 1) && !(TNF[u - 1] == 1 ||
                    IL6[u - 1] == 1 ||
                    IFNG == 1
                    ))
                {
                    if (p7 <= noiseLimit) { PPARg.Add(0); } else { PPARg.Add(1); }
                }
                else
                {
                    if (p7 <= noiseLimit) { PPARg.Add(1); } else { PPARg.Add(0); }
                }


                //Logic rule for CTGF
                if ((ROS[u - 1] == 1) ||
                    (CTGF[u - 1] == 1 && CTGFR == 1))
                {
                    if (p8 <= noiseLimit) { CTGF.Add(0); } else { CTGF.Add(1); }
                }
                else
                {
                    if (p8 <= noiseLimit) { CTGF.Add(1); } else { CTGF.Add(0); }
                }


                //Logic rule for Akt
                if ((Insulin == 1) && Ceramide[u - 1] == 0)
                {
                    if (p9 <= noiseLimit) { Akt.Add(0); } else { Akt.Add(1); }
                }
                else
                {
                    if (p9 <= noiseLimit) { Akt.Add(1); } else { Akt.Add(0); }
                }

                //Logic rule for Resistin
                if (TNF[u - 1] == 1 && AP1[u - 1] == 1 && PPARg[u - 1] == 0)
                {
                    if (p10 <= noiseLimit) { Resistin.Add(0); } else { Resistin.Add(1); }
                }
                else
                {
                    if (p10 <= noiseLimit) { Resistin.Add(1); } else { Resistin.Add(0); }
                }

                //Logic rule for Adiponectin
                if (PPARg[u - 1] == 1 && TNF[u - 1] == 0)
                {
                    if (p11 <= noiseLimit) { Adiponectin.Add(0); } else { Adiponectin.Add(1); }
                }
                else
                {
                    if (p11 <= noiseLimit) { Adiponectin.Add(1); } else { Adiponectin.Add(0); }
                }

                //Logic rule for Ceramide
                if ((sFFA[u - 1] == 1 || ROS[u - 1] == 1 || TNF[u - 1] == 1 || CeramideExt == 1) && Adiponectin[u - 1] == 0)
                {
                    if (p12 <= noiseLimit) { Ceramide.Add(0); } else { Ceramide.Add(0); }
                }
                else
                {
                    if (p12 <= noiseLimit) { Ceramide.Add(0); } else { Ceramide.Add(0); }
                }

                //Logic rule for ROS
                if ((PPARg[u - 1] == 0 || MechanicalPressure == 1) && Ceramide[u - 1] == 1)
                {
                    if (p13 <= noiseLimit) { ROS.Add(0); } else { ROS.Add(1); }
                }
                else
                {
                    if (p13 <= noiseLimit) { ROS.Add(1); } else { ROS.Add(0); }
                }

                //Logic rule for sFFA
                if (Akt[u - 1] == 0 && Ceramide[u - 1] == 1)
                {
                    if (p14 <= noiseLimit) { sFFA.Add(0); } else { sFFA.Add(1); }
                }
                else
                {
                    if (p14 <= noiseLimit) { sFFA.Add(1); } else { sFFA.Add(0); }
                }

                //Logic rule for GLUT4
                if (Akt[u - 1] == 1 && Ceramide[u - 1] == 0)
                {
                    if (p15 <= noiseLimit) { GLUT4.Add(0); } else { GLUT4.Add(1); }
                }
                else
                {
                    if (p14 <= noiseLimit) { GLUT4.Add(1); } else { GLUT4.Add(0); }
                }
            }


            //Subprocess 3: Saving data
            var node1 = IL1b.Select(x => Convert.ToString(x)).ToList();
            var node2 = IL6.Select(x => Convert.ToString(x)).ToList();
            var node3 = IL8.Select(x => Convert.ToString(x)).ToList();
            var node4 = TNF.Select(x => Convert.ToString(x)).ToList();
            var node5 = IL10.Select(x => Convert.ToString(x)).ToList();
            var node6 = IL4.Select(x => Convert.ToString(x)).ToList();
            var node7 = AP1.Select(x => Convert.ToString(x)).ToList();
            var node8 = PPARg.Select(x => Convert.ToString(x)).ToList();
            var node9 = CTGF.Select(x => Convert.ToString(x)).ToList();
            var node10 = Akt.Select(x => Convert.ToString(x)).ToList();
            var node11 = Resistin.Select(x => Convert.ToString(x)).ToList();
            var node12 = Adiponectin.Select(x => Convert.ToString(x)).ToList();
            var node13 = Ceramide.Select(x => Convert.ToString(x)).ToList();
            var node14 = ROS.Select(x => Convert.ToString(x)).ToList();
            var node15 = sFFA.Select(x => Convert.ToString(x)).ToList();
            var node16 = GLUT4.Select(x => Convert.ToString(x)).ToList();


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
                                    node15[w] + "  " +
                                    node16[w]);

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
        int xTNF_xCTGF_xGLUT4 = 0;
        int xTNF_xCTGF_GLUT4 = 0;
        int xTNF_CTGF_xGLUT4 = 0;
        int TNF_xCTGF_xGLUT4 = 0;
        int xTNF_CTGF_GLUT4 = 0;
        int TNF_CTGF_xGLUT4 = 0;
        int TNF_xCTGF_GLUT4 = 0;
        int TNF_CTGF_GLUT4 = 0;
        int totalTNF = 0;
        int totalCTGF = 0;
        int totalGLUT4 = 0;
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
        outcomeModel.Write("      Nodes ID: {0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15} Initiated in  : {" + intialConditions[hn] + "}: " + CellTags[hn] + "\n");

        foreach (var x in q)
        {
            Console.WriteLine("Fixed point in: {" + x.Value + "}  Size of its basin of attraction: {" + x.Count + "}");
            outcomeModel.Write("Fixed point in: {" + x.Value + "}  Size of its basin of attraction: {" + x.Count + "}\n");
            //Labelling rules
            string ste = x.Value.Replace("  ", "");
            char[] Repase = ste.ToCharArray();
            Console.WriteLine(Repase.Length);

            if (Repase[3] == '0' && Repase[8] == '0' && Repase[15] == '0')
            {
                xTNF_xCTGF_xGLUT4 = xTNF_xCTGF_xGLUT4 + x.Count;
            }
            if (Repase[3] == '0' && Repase[8] == '0' && Repase[15] == '1')
            {
                xTNF_xCTGF_GLUT4 = xTNF_xCTGF_GLUT4 + x.Count;
            }
            if (Repase[3] == '0' && Repase[8] == '1' && Repase[15] == '0')
            {
                xTNF_CTGF_xGLUT4 = xTNF_CTGF_xGLUT4 + x.Count;
            }
            if (Repase[3] == '1' && Repase[8] == '0' && Repase[15] == '0')
            {
                TNF_xCTGF_xGLUT4 = TNF_xCTGF_xGLUT4 + x.Count;
            }
            if (Repase[3] == '0' && Repase[8] == '1' && Repase[15] == '1')
            {
                xTNF_CTGF_GLUT4 = xTNF_CTGF_GLUT4 + x.Count;
            }
            if (Repase[3] == '0' && Repase[8] == '1' && Repase[15] == '1')
            {
                xTNF_CTGF_GLUT4 = xTNF_CTGF_GLUT4 + x.Count;
            }
            if (Repase[3] == '1' && Repase[8] == '0' && Repase[15] == '1')
            {
                TNF_xCTGF_GLUT4 = TNF_xCTGF_GLUT4 + x.Count;
            }
            if (Repase[3] == '1' && Repase[8] == '1' && Repase[15] == '0')
            {
                TNF_CTGF_xGLUT4 = TNF_CTGF_xGLUT4 + x.Count;
            }
            if (Repase[3] == '1' && Repase[8] == '1' && Repase[15] == '1')
            {
                TNF_CTGF_GLUT4 = TNF_CTGF_GLUT4 + x.Count;
            }

            if (Repase[3] == '1')
            {
                totalTNF = totalTNF + x.Count;
            }

            if (Repase[8] == '1')
            {
                totalCTGF = totalCTGF + x.Count;
            }

            if (Repase[15] == '1')
            {
                totalGLUT4 = totalGLUT4 + x.Count;
            }



        }
        outcomeModel.Write("\n");


        //write
        totalEnd = xTNF_xCTGF_xGLUT4 + xTNF_xCTGF_GLUT4 + xTNF_CTGF_xGLUT4 + TNF_xCTGF_xGLUT4 + xTNF_CTGF_GLUT4 + TNF_CTGF_xGLUT4 + TNF_xCTGF_GLUT4 + TNF_CTGF_GLUT4;

        outcomeModel.Write("\n");
        outcomeModel.Write("\n");
        outcomeModel.WriteLine("Total count: " + totalEnd);

        outcomeModel.WriteLine("TNF(-) CTGF(-) GLUT4(-): " + xTNF_xCTGF_xGLUT4);
        outcomeModel.WriteLine("TNF(-) CTGF(-) GLUT4(+): " + xTNF_xCTGF_GLUT4);
        outcomeModel.WriteLine("TNF(-) CTGF(+) GLUT4(-): " + xTNF_CTGF_xGLUT4);
        outcomeModel.WriteLine("TNF(-) CTGF(+) GLUT4(+): " + xTNF_CTGF_GLUT4);
        outcomeModel.WriteLine("TNF(+) CTGF(-) GLUT4(-): " + TNF_xCTGF_xGLUT4);
        outcomeModel.WriteLine("TNF(+) CTGF(-) GLUT4(+): " + TNF_xCTGF_GLUT4);
        outcomeModel.WriteLine("TNF(+) CTGF(+) GLUT4(-): " + TNF_CTGF_xGLUT4);
        outcomeModel.WriteLine("TNF(+) CTGF(+) GLUT4(+): " + TNF_CTGF_GLUT4);

        double resultxTNF_xCTGF_xGLUT4 = 0;
        double resultxTNF_xCTGF_GLUT4 = 0;
        double resultxTNF_CTGF_xGLUT4 = 0;
        double resultTNF_xCTGF_xGLUT4 = 0;
        double resultxTNF_CTGF_GLUT4 = 0;
        double resultTNF_xCTGF_GLUT4 = 0;
        double resultTNF_CTGF_xGLUT4 = 0;
        double resultTNF_CTGF_GLUT4 = 0;



        resultxTNF_xCTGF_xGLUT4 = (double)xTNF_xCTGF_xGLUT4 / (double)totalEnd;
        resultxTNF_xCTGF_GLUT4 = (double)xTNF_xCTGF_GLUT4 / (double)totalEnd;
        resultxTNF_CTGF_xGLUT4 = (double)xTNF_CTGF_xGLUT4 / (double)totalEnd;
        resultTNF_xCTGF_xGLUT4 = (double)TNF_xCTGF_xGLUT4 / (double)totalEnd;
        resultxTNF_CTGF_GLUT4 = (double)xTNF_CTGF_GLUT4 / (double)totalEnd;
        resultTNF_xCTGF_GLUT4 = (double)TNF_xCTGF_GLUT4 / (double)totalEnd;
        resultTNF_CTGF_xGLUT4 = (double)TNF_CTGF_xGLUT4 / (double)totalEnd;
        resultTNF_CTGF_GLUT4 = (double)TNF_CTGF_GLUT4 / (double)totalEnd;

        outcomeModel.Write("\n");


        onlyDataModel.WriteLine(resultxTNF_xCTGF_xGLUT4);
        onlyDataModel.WriteLine(resultxTNF_xCTGF_GLUT4);
        onlyDataModel.WriteLine(resultxTNF_CTGF_xGLUT4);
        onlyDataModel.WriteLine(resultxTNF_CTGF_GLUT4);
        onlyDataModel.WriteLine(resultTNF_xCTGF_xGLUT4);
        onlyDataModel.WriteLine(resultTNF_xCTGF_GLUT4);
        onlyDataModel.WriteLine(resultTNF_CTGF_xGLUT4);
        onlyDataModel.WriteLine(resultTNF_CTGF_GLUT4);

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
    outcomeModel.Write(" IL1b: {0}\n IL6: {1}\n IL8: {2}\n TNF: {3}\n IL10: {4}\n");
    outcomeModel.Write(" IL4: {5}\n AP1: {6}\n PPARg: {7}\n CTGF: {8}\n Akt: {9}\n");
    outcomeModel.Write(" Resistin: {10}\n Adiponectin: {11}\n Ceramide: {12}\n ROS: {13}\n sFFA: {14}\n GLUT4: {15}\n");

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
