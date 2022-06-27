// See https://aka.ms/new-console-template for more information
for (int k = 0; k <= 10; k++)
{
    string fileName = "CD4_T_cells";
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
    List<string> fixedPoints8 = new List<string>();
    List<string> fixedPoints9 = new List<string>();

    StreamWriter outcomeModel = new StreamWriter(filepath);
    StreamWriter onlyDataModel = new StreamWriter(filepath1);

    List<string>[] fixe = { fixedPoints, fixedPoints1, fixedPoints2, fixedPoints3, fixedPoints4, fixedPoints5, fixedPoints6, fixedPoints7, fixedPoints8, fixedPoints9 };

    int maxN = 10;
    int maxN1 = 10000;
    int maxNtotal = (int)Math.Pow(2, maxN);
    int iterationN = 30;
    double noiseLimit = 0.08;
    string[] CellTags = { "Th0", "Th1", "Th2", "Th9", "Th17", "Th1R", "Th2R", "iTreg", "Tr1", "Th3" };

    string a1 = "0  0  0  0  0  0  0  0  0  0  0  0";
    string a2 = "1  1  0  0  0  0  0  0  0  0  0  0";
    string a3 = "0  0  1  1  1  0  0  0  0  0  0  0";
    string a4 = "0  0  0  0  0  0  0  0  1  1  1  1";
    string a5 = "0  0  0  0  0  1  1  0  1  0  1  0";
    string a6 = "1  1  0  0  0  0  0  0  1  0  1  0";
    string a7 = "0  0  1  0  1  0  0  0  0  1  0  0";
    string a8 = "0  0  0  1  0  0  0  1  1  0  0  0";
    string a9 = "0  0  0  0  0  0  0  0  0  1  0  0";
    string a10 = "0  0  0  0  0  0  0  0  1  1  1  0";

    intialConditions.Add(a1.Replace("  ", ""));
    intialConditions.Add(a2.Replace("  ", ""));
    intialConditions.Add(a3.Replace("  ", ""));
    intialConditions.Add(a4.Replace("  ", ""));
    intialConditions.Add(a5.Replace("  ", ""));
    intialConditions.Add(a6.Replace("  ", ""));
    intialConditions.Add(a7.Replace("  ", ""));
    intialConditions.Add(a8.Replace("  ", ""));
    intialConditions.Add(a9.Replace("  ", ""));
    intialConditions.Add(a10.Replace("  ", ""));

    #endregion
    //*****************

    //Step 3: Setting input for cases of study
    int IFNGe = 0;
    int IL12e = 0;
    int IL2e = 0;
    int IL4e = 01;
    int TGFBe = 0;
    int IL10e = 01;
    int IL21e = 0;
    int IL27e = 0;
    int INSULIN = 01;
    int Ceramide = 01;


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
            List<int> TBET = new List<int>();
            List<int> IFNG = new List<int>();
            List<int> GATA3 = new List<int>();
            List<int> IL2 = new List<int>();
            List<int> IL4 = new List<int>();
            List<int> RORGT = new List<int>();
            List<int> IL6 = new List<int>();
            List<int> FOXP3 = new List<int>();
            List<int> TGFB = new List<int>();
            List<int> IL10 = new List<int>();
            List<int> PU1 = new List<int>();
            List<int> IL9 = new List<int>();

            TBET.Add(xS[0]);
            IFNG.Add(xS[1]);
            GATA3.Add(xS[2]);
            IL2.Add(xS[3]);
            IL4.Add(xS[4]);
            RORGT.Add(xS[5]);
            IL6.Add(xS[6]);
            FOXP3.Add(xS[7]);
            TGFB.Add(xS[8]);
            IL10.Add(xS[9]);
            PU1.Add(xS[10]);
            IL9.Add(xS[11]);


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

                //Logic rule for TBET
                if (((IFNG[u - 1] == 1 || (
                    IL12e == 1 && !(
                    IL6[u - 1] == 1 ||
                    IL4[u - 1] == 1 ||
                    IL10[u - 1] == 1))) ||
                    TBET[u - 1] == 1) && !(
                    IL4[u - 1] == 1 ||
                    GATA3[u - 1] == 1 ||
                    IL6[u - 1] == 1))
                {
                    if (p0 <= noiseLimit) { TBET.Add(0); } else { TBET.Add(1); }
                }
                else
                {
                    if (p0 <= noiseLimit) { TBET.Add(1); } else { TBET.Add(0); }

                }


                //Logic rule for IFNG
                if ((IFNGe == 1 || ((
                    IFNG[u - 1] == 1 ||
                    TBET[u - 1] == 1) && !(
                    GATA3[u - 1] == 1 ||
                    TGFB[u - 1] == 1))) && !(
                    IL6[u - 1] == 1 ||
                    IL4[u - 1] == 1 ||
                    IL10[u - 1] == 1 || Ceramide == 1))
                {
                    if (p1 <= noiseLimit) { IFNG.Add(0); } else { IFNG.Add(1); }

                }
                else
                {
                    if (p1 <= noiseLimit) { IFNG.Add(1); } else { IFNG.Add(0); }

                }


                //Logic rule for GATA3
                if (((IL2[u - 1] == 1 &&
                    IL4[u - 1] == 1) ||
                    GATA3[u - 1] == 1) && !(
                    TBET[u - 1] == 1 ||
                    TGFB[u - 1] == 1 ||
                    IL6[u - 1] == 1 ||
                    IFNG[u - 1] == 1 ||
                    PU1[u - 1] == 1))
                {
                    if (p2 <= noiseLimit) { GATA3.Add(0); } else { GATA3.Add(1); }
                }
                else
                {
                    if (p2 <= noiseLimit) { GATA3.Add(1); } else { GATA3.Add(0); }

                }


                //Logic rule for IL2
                if ((IL2e == 1 || (
                    IL2[u - 1] == 1 &&
                    FOXP3[u - 1] == 0)) && !(
                    IFNG[u - 1] == 1 ||
                    IL6[u - 1] == 1 || (
                    IL10[u - 1] == 1 &&
                    FOXP3[u - 1] == 0)))
                {
                    if (p3 <= noiseLimit) { IL2.Add(0); } else { IL2.Add(1); }
                }
                else
                {
                    if (p3 <= noiseLimit) { IL2.Add(1); } else { IL2.Add(0); }
                }


                //Logic rule for IL4
                if ((IL4e == 1 || (
                    GATA3[u - 1] == 1 && (
                    IL2[u - 1] == 1 ||
                    IL4[u - 1] == 1) &&
                    TBET[u - 1] == 0)) && !(
                    IFNG[u - 1] == 1 ||
                    IL6[u - 1] == 1 ||
                    PU1[u - 1] == 1))
                {
                    if (p4 <= noiseLimit) { IL4.Add(0); } else { IL4.Add(1); }
                }
                else
                {
                    if (p4 <= noiseLimit) { IL4.Add(1); } else { IL4.Add(0); }
                }


                //Logic rule for RORGT
                if ((IL6[u - 1] == 1 &&
                    TGFB[u - 1] == 1) && !(
                    TBET[u - 1] == 1 ||
                    FOXP3[u - 1] == 1 ||
                    GATA3[u - 1] == 1))
                {
                    if (p5 <= noiseLimit) { RORGT.Add(0); } else { RORGT.Add(1); }
                }
                else
                {
                    if (p5 <= noiseLimit) { RORGT.Add(1); } else { RORGT.Add(0); }
                }


                //Logic rule for IL6
                if ((IL21e == 1 ||
                    IL6[u - 1] == 1 ||
                    RORGT[u - 1] == 1) && !(
                    IFNG[u - 1] == 1 ||
                    IL4[u - 1] == 1 ||
                    IL10[u - 1] == 1 ||
                    IL2[u - 1] == 1))
                {
                    if (p6 <= noiseLimit) { IL6.Add(0); } else { IL6.Add(1); }
                }
                else
                {
                    if (p6 <= noiseLimit) { IL6.Add(1); } else { IL6.Add(0); }
                }


                //Logic rule for FOXP3
                if ((IL2[u - 1] == 1 && (Ceramide == 1 ||
                    TGFB[u - 1] == 1 ||
                    FOXP3[u - 1] == 1)) && !(
                    IL6[u - 1] == 1 ||
                    RORGT[u - 1] == 1 ||
                    PU1[u - 1] == 1))
                {
                    if (p7 <= noiseLimit) { FOXP3.Add(0); } else { FOXP3.Add(1); }
                }
                else
                {
                    if (p7 <= noiseLimit) { FOXP3.Add(1); } else { FOXP3.Add(0); }
                }


                //Logic rule for TGFB
                if (TGFBe == 1 || ((
                    TGFB[u - 1] == 1 ||
                    FOXP3[u - 1] == 1) &&
                    IL6[u - 1] == 0))
                {
                    if (p8 <= noiseLimit) { TGFB.Add(0); } else { TGFB.Add(1); }
                }
                else
                {
                    if (p8 <= noiseLimit) { TGFB.Add(1); } else { TGFB.Add(0); }
                }


                //Logic rule for IL10
                if (IL10e == 1 || (
                    IL10[u - 1] == 1 && (
                    IFNG[u - 1] == 1 ||
                    IL6[u - 1] == 1 ||
                    TGFB[u - 1] == 1 ||
                    GATA3[u - 1] == 1 ||
                    IL27e == 1) &&
                    !(INSULIN == 1 || Ceramide == 1)))
                {
                    if (p9 <= noiseLimit) { IL10.Add(0); } else { IL10.Add(1); }
                }
                else
                {
                    if (p9 <= noiseLimit) { IL10.Add(1); } else { IL10.Add(0); }
                }

                //Logic rule for PU1
                if ((IL6[u - 1] == 1 ||
                    TGFB[u - 1] == 1 ||
                    PU1[u - 1] == 1) &&
                    IL2[u - 1] == 0)
                {
                    if (p10 <= noiseLimit) { PU1.Add(0); } else { PU1.Add(1); }
                }
                else
                {
                    if (p10 <= noiseLimit) { PU1.Add(1); } else { PU1.Add(0); }
                }

                //Logic rule for IL9
                if (((IL4[u - 1] == 1 ||
                    IL4e == 1) &&
                    PU1[u - 1] == 1) &&
                    TBET[u - 1] == 0)
                {
                    if (p11 <= noiseLimit) { IL9.Add(0); } else { IL9.Add(1); }
                }
                else
                {
                    if (p11 <= noiseLimit) { IL9.Add(1); } else { IL9.Add(0); }
                }

            }


            //Subprocess 3: Saving data
            var node1 = TBET.Select(x => Convert.ToString(x)).ToList();
            var node2 = IFNG.Select(x => Convert.ToString(x)).ToList();
            var node3 = GATA3.Select(x => Convert.ToString(x)).ToList();
            var node4 = IL2.Select(x => Convert.ToString(x)).ToList();
            var node5 = IL4.Select(x => Convert.ToString(x)).ToList();
            var node6 = RORGT.Select(x => Convert.ToString(x)).ToList();
            var node7 = IL6.Select(x => Convert.ToString(x)).ToList();
            var node8 = FOXP3.Select(x => Convert.ToString(x)).ToList();
            var node9 = TGFB.Select(x => Convert.ToString(x)).ToList();
            var node10 = IL10.Select(x => Convert.ToString(x)).ToList();
            var node11 = PU1.Select(x => Convert.ToString(x)).ToList();
            var node12 = IL9.Select(x => Convert.ToString(x)).ToList();


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
                          node12[w]);

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
        int Th0 = 0;
        int Th1 = 0;
        int Th2 = 0;
        int Th9 = 0;
        int Th17 = 0;
        int Th1R = 0;
        int Th2R = 0;
        int iTreg = 0;
        int Tr1 = 0;
        int Th3 = 0;
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
        outcomeModel.Write("      Nodes ID: {0  1  2  3  4  5  6  7  8  9  10 11} Initiated in  : {" + intialConditions[hn] + "}: " + CellTags[hn] + "\n");

        foreach (var x in q)
        {
            Console.WriteLine("Fixed point in: {" + x.Value + "}  Size of its basin of attraction: {" + x.Count + "}");
            outcomeModel.Write("Fixed point in: {" + x.Value + "}  Size of its basin of attraction: {" + x.Count + "}\n");
            //Labelling rules
            string ste = x.Value.Replace("  ", "");
            char[] Repase = ste.ToCharArray();
            Console.WriteLine(Repase.Length);

            if (!(Repase[0] == '1' ||
                Repase[2] == '1' ||
                Repase[5] == '1' ||
                Repase[7] == '1' ||
                Repase[9] == '1' ||
                Repase[8] == '1' ||
                Repase[11] == '1'))
            {
                Th0 = Th0 + x.Count;
            }
            if ((Repase[0] == '1' &&
                    Repase[1] == '1') && !(
                    Repase[9] == '1' ||
                    Repase[8] == '1' ||
                    Repase[7] == '1'))
            {
                Th1 = Th1 + x.Count;
            }
            if ((Repase[2] == '1' &&
                    Repase[4] == '1') && !(
                    Repase[9] == '1' ||
                    Repase[8] == '1' ||
                    Repase[7] == '1'))
            {
                Th2 = Th2 + x.Count;
            }
            if ((Repase[10] == '1' &&
                    Repase[11] == '1') && !(
                    Repase[4] == '1' ||
                    Repase[1] == '1' ||
                    Repase[0] == '1' ||
                    Repase[2] == '1' ||
                    Repase[5] == '1' ||
                    Repase[7] == '1'))
            {
                Th9 = Th9 + x.Count;
            }
            if (Repase[5] == '1' &&
                    Repase[6] == '1' &&
                    Repase[9] == '0')
            {
                Th17 = Th17 + x.Count;
            }
            if (Repase[0] == '1' && (
                    Repase[9] == '1' ||
                    Repase[8] == '1' ||
                    Repase[7] == '1'))
            {
                Th1R = Th1R + x.Count;
            }
            if (Repase[2] == '1' && (
                    Repase[9] == '1' ||
                    Repase[8] == '1' ||
                    Repase[7] == '1'))
            {
                Th2R = Th2R + x.Count;
            }
            if (Repase[7] == '1' &&
                    Repase[8] == '1' && !(
                    Repase[0] == '1' ||
                    Repase[2] == '1' ||
                    Repase[5] == '1'))
            {
                iTreg = iTreg + x.Count;
            }
            if (Repase[9] == '1' && !(
                    Repase[8] == '1' ||
                    Repase[2] == '1' ||
                    Repase[7] == '1' ||
                    Repase[5] == '1'))
            {
                Tr1 = Tr1 + x.Count;
            }
            if (Repase[8] == '1' && !(
                    Repase[0] == '1' ||
                    Repase[2] == '1' ||
                    Repase[7] == '1' ||
                    Repase[5] == '1'))
            {
                Th3 = Th3 + x.Count;
            }


        }
        outcomeModel.Write("\n");


        //write

        outcomeModel.Write("\n");
        outcomeModel.Write("\n");
        totalEnd = Th0 + Th1 + Th2 + Th9 + Th17 + Th1R + Th2R + iTreg + Tr1 + Th3;
        outcomeModel.WriteLine("Total count: " + totalEnd);
        double resultTh0 = 0;
        double resultTh1 = 0;
        double resultTh2 = 0;
        double resultTh9 = 0;
        double resultTh17 = 0;
        double resultTh1R = 0;
        double resultTh2R = 0;
        double resultiTreg = 0;
        double resultTr1 = 0;
        double resultTh3 = 0;

        resultTh0 = (double)Th0 / (double)totalEnd;
        resultTh1 = (double)Th1 / (double)totalEnd;
        resultTh2 = (double)Th2 / (double)totalEnd;
        resultTh9 = (double)Th9 / (double)totalEnd;
        resultTh17 = (double)Th17 / (double)totalEnd;
        resultTh1R = (double)Th1R / (double)totalEnd;
        resultTh2R = (double)Th2R / (double)totalEnd;
        resultiTreg = (double)iTreg / (double)totalEnd;
        resultTr1 = (double)Tr1 / (double)totalEnd;
        resultTh3 = (double)Th3 / (double)totalEnd;

        outcomeModel.Write("\n");

        outcomeModel.WriteLine("Th0: " + Th0);
        outcomeModel.WriteLine("Th1: " + Th1);
        outcomeModel.WriteLine("Th2: " + Th2);
        outcomeModel.WriteLine("Th9: " + Th9);
        outcomeModel.WriteLine("Th17: " + Th17);
        outcomeModel.WriteLine("Th1R: " + Th2R);
        outcomeModel.WriteLine("Th2R: " + Th1R);
        outcomeModel.WriteLine("iTreg: " + iTreg);
        outcomeModel.WriteLine("Tr1: " + Tr1);
        outcomeModel.WriteLine("Th3: " + Th3);
        outcomeModel.Write("\n");

        outcomeModel.WriteLine(Th0);
        outcomeModel.WriteLine(Th1);
        outcomeModel.WriteLine(Th2);
        outcomeModel.WriteLine(Th9);
        outcomeModel.WriteLine(Th17);
        outcomeModel.WriteLine(Th2R);
        outcomeModel.WriteLine(Th1R);
        outcomeModel.WriteLine(iTreg);
        outcomeModel.WriteLine(Tr1);
        outcomeModel.WriteLine(Th3);


        onlyDataModel.WriteLine(resultTh0);
        onlyDataModel.WriteLine(resultTh1);
        onlyDataModel.WriteLine(resultTh2);
        onlyDataModel.WriteLine(resultTh9);
        onlyDataModel.WriteLine(resultTh17);
        onlyDataModel.WriteLine(resultTh1R);
        onlyDataModel.WriteLine(resultTh2R);
        onlyDataModel.WriteLine(resultiTreg);
        onlyDataModel.WriteLine(resultTr1);
        onlyDataModel.WriteLine(resultTh3);

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
    outcomeModel.Write(" TBET: {0}\n IFNG: {1}\n GATA3: {2}\n IL2: {3}\n IL4: {4}\n");
    outcomeModel.Write(" RORGT: {5}\n IL6: {6}\n FOXP3: {7}\n TGFB: {8}\n IL10: {9}\n");
    outcomeModel.Write(" PU1: {10}\n IL9: {11}\n");
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
