using System;
using System.Drawing;
using System.Windows.Forms;

namespace Project2
{
    public partial class Form1 : Form
    {
        public bool epipolar = false;
        public static void svdcmp(double[,] a, out double[] w, out double[,] v)
        {
            int m = a.GetLength(0);
            // number of columns in A
            int n = a.GetLength(1);
            /*if (m < n)
            {
                throw new ArgumentException("Number of rows in A must be greater or equal to number of columns");
            }*/

            w = new double[n];
            v = new double[n, n];


            int flag, i, its, j, jj, k, l = 0, nm = 0;
            double anorm, c, f, g, h, s, scale, x, y, z;

            double[] rv1 = new double[n];

            // householder reduction to bidiagonal form
            g = scale = anorm = 0.0;

            for (i = 0; i < n; i++)
            {
                l = i + 1;
                rv1[i] = scale * g;
                g = s = scale = 0;

                if (i < m)
                {
                    for (k = i; k < m; k++)
                    {
                        scale += System.Math.Abs(a[k, i]);
                    }

                    if (scale != 0.0)
                    {
                        for (k = i; k < m; k++)
                        {
                            a[k, i] /= scale;
                            s += a[k, i] * a[k, i];
                        }

                        f = a[i, i];
                        g = -Sign(System.Math.Sqrt(s), f);
                        h = f * g - s;
                        a[i, i] = f - g;

                        if (i != n - 1)
                        {
                            for (j = l; j < n; j++)
                            {
                                for (s = 0.0, k = i; k < m; k++)
                                {
                                    s += a[k, i] * a[k, j];
                                }

                                f = s / h;

                                for (k = i; k < m; k++)
                                {
                                    a[k, j] += f * a[k, i];
                                }
                            }
                        }

                        for (k = i; k < m; k++)
                        {
                            a[k, i] *= scale;
                        }
                    }
                }

                w[i] = scale * g;
                g = s = scale = 0.0;

                if ((i < m) && (i != n - 1))
                {
                    for (k = l; k < n; k++)
                    {
                        scale += System.Math.Abs(a[i, k]);
                    }

                    if (scale != 0.0)
                    {
                        for (k = l; k < n; k++)
                        {
                            a[i, k] /= scale;
                            s += a[i, k] * a[i, k];
                        }

                        f = a[i, l];
                        g = -Sign(System.Math.Sqrt(s), f);
                        h = f * g - s;
                        a[i, l] = f - g;

                        for (k = l; k < n; k++)
                        {
                            rv1[k] = a[i, k] / h;
                        }

                        if (i != m - 1)
                        {
                            for (j = l; j < m; j++)
                            {
                                for (s = 0.0, k = l; k < n; k++)
                                {
                                    s += a[j, k] * a[i, k];
                                }
                                for (k = l; k < n; k++)
                                {
                                    a[j, k] += s * rv1[k];
                                }
                            }
                        }

                        for (k = l; k < n; k++)
                        {
                            a[i, k] *= scale;
                        }
                    }
                }
                anorm = System.Math.Max(anorm, (System.Math.Abs(w[i]) + System.Math.Abs(rv1[i])));
            }

            // accumulation of right-hand transformations
            for (i = n - 1; i >= 0; i--)
            {
                if (i < n - 1)
                {
                    if (g != 0.0)
                    {
                        for (j = l; j < n; j++)
                        {
                            v[j, i] = (a[i, j] / a[i, l]) / g;
                        }

                        for (j = l; j < n; j++)
                        {
                            for (s = 0, k = l; k < n; k++)
                            {
                                s += a[i, k] * v[k, j];
                            }
                            for (k = l; k < n; k++)
                            {
                                v[k, j] += s * v[k, i];
                            }
                        }
                    }
                    for (j = l; j < n; j++)
                    {
                        v[i, j] = v[j, i] = 0;
                    }
                }
                v[i, i] = 1;
                g = rv1[i];
                l = i;
            }

            // accumulation of left-hand transformations
            for (i = n - 1; i >= 0; i--)
            {
                l = i + 1;
                g = w[i];

                if (i < n - 1)
                {
                    for (j = l; j < n; j++)
                    {
                        a[i, j] = 0.0;
                    }
                }

                if (g != 0)
                {
                    g = 1.0 / g;

                    if (i != n - 1)
                    {
                        for (j = l; j < n; j++)
                        {
                            for (s = 0, k = l; k < m; k++)
                            {
                                s += a[k, i] * a[k, j];
                            }

                            f = (s / a[i, i]) * g;

                            for (k = i; k < m; k++)
                            {
                                a[k, j] += f * a[k, i];
                            }
                        }
                    }

                    for (j = i; j < m; j++)
                    {
                        a[j, i] *= g;
                    }
                }
                else
                {
                    for (j = i; j < m; j++)
                    {
                        a[j, i] = 0;
                    }
                }
                ++a[i, i];
            }

            // diagonalization of the bidiagonal form: Loop over singular values
            // and over allowed iterations
            for (k = n - 1; k >= 0; k--)
            {
                for (its = 1; its <= 30; its++)
                {
                    flag = 1;

                    for (l = k; l >= 0; l--)
                    {
                        // test for splitting
                        nm = l - 1;

                        if (System.Math.Abs(rv1[l]) + anorm == anorm)
                        {
                            flag = 0;
                            break;
                        }

                        if (System.Math.Abs(w[nm]) + anorm == anorm)
                            break;
                    }

                    if (flag != 0)
                    {
                        c = 0.0;
                        s = 1.0;
                        for (i = l; i <= k; i++)
                        {
                            f = s * rv1[i];

                            if (System.Math.Abs(f) + anorm != anorm)
                            {
                                g = w[i];
                                h = Pythag(f, g);
                                w[i] = h;
                                h = 1.0 / h;
                                c = g * h;
                                s = -f * h;

                                for (j = 0; j < m; j++)
                                {
                                    y = a[j, nm];
                                    z = a[j, i];
                                    a[j, nm] = y * c + z * s;
                                    a[j, i] = z * c - y * s;
                                }
                            }
                        }
                    }

                    z = w[k];

                    if (l == k)
                    {
                        // convergence
                        if (z < 0.0)
                        {
                            // singular value is made nonnegative
                            w[k] = -z;

                            for (j = 0; j < n; j++)
                            {
                                v[j, k] = -v[j, k];
                            }
                        }
                        break;
                    }

                    if (its == 30)
                    {
                        throw new ApplicationException("No convergence in 30 svdcmp iterations");
                    }

                    // shift from bottom 2-by-2 minor
                    x = w[l];
                    nm = k - 1;
                    y = w[nm];
                    g = rv1[nm];
                    h = rv1[k];
                    f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
                    g = Pythag(f, 1.0);
                    f = ((x - z) * (x + z) + h * ((y / (f + Sign(g, f))) - h)) / x;

                    // next QR transformation
                    c = s = 1.0;

                    for (j = l; j <= nm; j++)
                    {
                        i = j + 1;
                        g = rv1[i];
                        y = w[i];
                        h = s * g;
                        g = c * g;
                        z = Pythag(f, h);
                        rv1[j] = z;
                        c = f / z;
                        s = h / z;
                        f = x * c + g * s;
                        g = g * c - x * s;
                        h = y * s;
                        y *= c;

                        for (jj = 0; jj < n; jj++)
                        {
                            x = v[jj, j];
                            z = v[jj, i];
                            v[jj, j] = x * c + z * s;
                            v[jj, i] = z * c - x * s;
                        }

                        z = Pythag(f, h);
                        w[j] = z;

                        if (z != 0)
                        {
                            z = 1.0 / z;
                            c = f * z;
                            s = h * z;
                        }

                        f = c * g + s * y;
                        x = c * y - s * g;

                        for (jj = 0; jj < m; jj++)
                        {
                            y = a[jj, j];
                            z = a[jj, i];
                            a[jj, j] = y * c + z * s;
                            a[jj, i] = z * c - y * s;
                        }
                    }

                    rv1[l] = 0.0;
                    rv1[k] = f;
                    w[k] = x;
                }
            }
        }

        private static double Sign(double a, double b)
        {
            return (b >= 0.0) ? System.Math.Abs(a) : -System.Math.Abs(a);
        }

        private static double Pythag(double a, double b)
        {
            double at = System.Math.Abs(a), bt = System.Math.Abs(b), ct, result;

            if (at > bt)
            {
                ct = bt / at;
                result = at * System.Math.Sqrt(1.0 + ct * ct);
            }
            else if (bt > 0.0)
            {
                ct = at / bt;
                result = bt * System.Math.Sqrt(1.0 + ct * ct);
            }
            else
            {
                result = 0.0;
            }

            return result;
        }

        public Graphics g;
        Pen p;
        double[,] F = new double[3, 3];
        Image defaultImage1, gImage1, defaultImage2, gImage2;
        public Form1()
        {
            InitializeComponent();
            textBox1.Visible = false;
            g = this.CreateGraphics();
            p = new Pen(Color.LightGreen, 1);
        }

        private void label1_Click(object sender, EventArgs e)
        {

        }

        

        private void button2_Click(object sender, EventArgs e)
        {
            OpenFileDialog open = new OpenFileDialog();
            // image filters  
            open.Filter = "Image Files(*.jpg; *.jpeg; *.gif; *.bmp)|*.jpg; *.jpeg; *.gif; *.bmp";
            if (open.ShowDialog() == DialogResult.OK)
            {
                // display image in picture box  
                pictureBox1.Image = new Bitmap(open.FileName);
                defaultImage1 = pictureBox1.Image;
                // image file path  
            }
        }

        private void button3_Click(object sender, EventArgs e)
        {
            OpenFileDialog open = new OpenFileDialog();
            // image filters  
            open.Filter = "Image Files(*.jpg; *.jpeg; *.gif; *.bmp)|*.jpg; *.jpeg; *.gif; *.bmp";
            if (open.ShowDialog() == DialogResult.OK)
            {
                // display image in picture box  
                pictureBox2.Image = new Bitmap(open.FileName);
                defaultImage2 = pictureBox2.Image;
                // image file path  
            }
        }
        public int[,] P = new int[20,2];
        public int[,] PD = new int[20, 2];
        public int i = 0;
        public double[] X = new double[3];
        public double[] Xd = new double[3];

        private void pictureBox1_Click(object sender, EventArgs e)
        {
            if(i%2 == 0)
            {
                MouseEventArgs me = (MouseEventArgs)e;
                Point coordinates = me.Location;

                if(epipolar == false)
                    listBox1.Items.Insert(i, "  " +  (i/2 + 1) + " => " + coordinates.ToString() + ",");
                
                P[i / 2, 0] = coordinates.X;
                P[i / 2, 1] = coordinates.Y;

                Graphics dots = pictureBox1.CreateGraphics();
                dots.FillRectangle(Brushes.LightGreen, coordinates.X - 3, coordinates.Y - 3, 6.0f, 6.0f);
                dots.DrawRectangle(Pens.Red, coordinates.X - 3, coordinates.Y - 3, 6.0f, 6.0f);

                /*if (epipolar)
                {
                    *//*
                        00 01 02    0
                        10 11 12    1
                        20 21 22    2

                   *//*
                    X[0] = P[i / 2, 0];
                    X[1] = P[i / 2, 1];
                    X[2] = 1;

                    double[] Eq = new double[3];

                    Eq[0] = ((F[0, 0] * X[0]) + (F[0, 1] * X[1]) + (F[0, 2] * X[2]));
                    Eq[1] = ((F[1, 0] * X[0]) + (F[1, 1] * X[1]) + (F[1, 2] * X[2]));
                    Eq[2] = ((F[2, 0] * X[0]) + (F[2, 1] * X[1]) + (F[2, 2] * X[2]));

                    int x1 = 0;
                    int x2 = 430;

                    double y1, y2;

                    y1 = (((-1) * Eq[2]) - (Eq[0] * x1)) / Eq[1];
                    y2 = (((-1) * Eq[2]) - (Eq[0] * x2)) / Eq[1];


                    Graphics lines = pictureBox2.CreateGraphics();
                    lines.DrawLine(Pens.LightGreen, x1, Convert.ToInt32(y1), x2, Convert.ToInt32(y2));
                }*/

                i = i+1;
            }

            if (epipolar)
            {
                MouseEventArgs me = (MouseEventArgs)e;
                Point coordinates = me.Location;

                if (epipolar == false)
                    listBox1.Items.Insert(i, "  " + (i / 2 + 1) + " => " + coordinates.ToString() + ",");

                P[i / 2, 0] = coordinates.X;
                P[i / 2, 1] = coordinates.Y;

                Graphics dots = pictureBox1.CreateGraphics();
                dots.FillRectangle(Brushes.LightGreen, coordinates.X - 3, coordinates.Y - 3, 6.0f, 6.0f);
                dots.DrawRectangle(Pens.Red, coordinates.X - 3, coordinates.Y - 3, 6.0f, 6.0f);

                
                
                    /*
                        00 01 02    0
                        10 11 12    1
                        20 21 22    2
                    F * P
                   */
                    X[0] = P[i / 2, 0];
                    X[1] = P[i / 2, 1];
                    X[2] = 1;

                    double[] Eq = new double[3];

                    Eq[0] = ((F[0, 0] * X[0]) + (F[0, 1] * X[1]) + (F[0, 2] * X[2]));
                    Eq[1] = ((F[1, 0] * X[0]) + (F[1, 1] * X[1]) + (F[1, 2] * X[2]));
                    Eq[2] = ((F[2, 0] * X[0]) + (F[2, 1] * X[1]) + (F[2, 2] * X[2]));

                    int x1 = 0;
                    int x2 = 430;

                    double y1, y2;

                    y1 = (((-1) * Eq[2]) - (Eq[0] * x1)) / Eq[1];
                    y2 = (((-1) * Eq[2]) - (Eq[0] * x2)) / Eq[1];


                    Graphics lines = pictureBox2.CreateGraphics();
                    lines.DrawLine(Pens.LightGreen, x1, Convert.ToInt32(y1), x2, Convert.ToInt32(y2));
                

                i = i + 1;
            }
        }


        private void pictureBox2_Click(object sender, EventArgs e)
        {
            if (i % 2 != 0)
            {
                MouseEventArgs me = (MouseEventArgs)e;
                Point coordinates = me.Location;
                Graphics dots = pictureBox2.CreateGraphics();
                dots.FillRectangle(Brushes.LightGreen, coordinates.X - 3 , coordinates.Y - 3, 6.0f, 6.0f);
                dots.DrawRectangle(Pens.Red, coordinates.X - 3, coordinates.Y - 3, 6.0f, 6.0f);

                if(epipolar == false)
                    listBox1.Items.Add("       " + coordinates.ToString());
     
                PD[i/2, 0] = coordinates.X;
                PD[i/2, 1] = coordinates.Y;

                /*if (epipolar)
                {
                    *//*
                       0    1   2       00 01 02    
                                        10 11 12    
                                        20 21 22    

                   *//*
                    X[0] = PD[i / 2, 0];
                    X[1] = PD[i / 2, 1];
                    X[2] = 1;

                    double[] Eq = new double[3];

                    Eq[0] = ((F[0, 0] * X[0]) + (F[1, 0] * X[1]) + (F[2, 0] * X[2]));
                    Eq[1] = ((F[0, 1] * X[0]) + (F[1, 1] * X[1]) + (F[2, 1] * X[2]));
                    Eq[2] = ((F[0, 2] * X[0]) + (F[1, 2] * X[1]) + (F[2, 2] * X[2]));

                    int x1 = 0;
                    int x2 = 430;

                    double y1, y2;

                    y1 = (((-1) * Eq[2]) - (Eq[0] * x1)) / Eq[1];
                    y2 = (((-1) * Eq[2]) - (Eq[0] * x2)) / Eq[1];


                    Graphics lines = pictureBox1.CreateGraphics();
                    lines.DrawLine(Pens.Yellow, x1, Convert.ToInt32(y1), x2, Convert.ToInt32(y2));
                }*/
                i = i + 1;
            }

            if (epipolar)
            {
                MouseEventArgs me = (MouseEventArgs)e;
                Point coordinates = me.Location;
                Graphics dots = pictureBox2.CreateGraphics();
                dots.FillRectangle(Brushes.LightGreen, coordinates.X - 3, coordinates.Y - 3, 6.0f, 6.0f);
                dots.DrawRectangle(Pens.Red, coordinates.X - 3, coordinates.Y - 3, 6.0f, 6.0f);

                if (epipolar == false)
                    listBox1.Items.Add("       " + coordinates.ToString());

                PD[i / 2, 0] = coordinates.X;
                PD[i / 2, 1] = coordinates.Y;

                
                    /*
                       0    1   2       00 01 02    
                                        10 11 12    
                                        20 21 22   
                    

                P[x,y]
                    F[]


                   */
                    X[0] = PD[i / 2, 0];
                    X[1] = PD[i / 2, 1];
                    X[2] = 1;

                    double[] Eq = new double[3];

                    Eq[0] = ((F[0, 0] * X[0]) + (F[1, 0] * X[1]) + (F[2, 0] * X[2]));
                    Eq[1] = ((F[0, 1] * X[0]) + (F[1, 1] * X[1]) + (F[2, 1] * X[2]));
                    Eq[2] = ((F[0, 2] * X[0]) + (F[1, 2] * X[1]) + (F[2, 2] * X[2]));

                    int x1 = 0;
                    int x2 = 430;

                    double y1, y2;

                    y1 = (((-1) * Eq[2]) - (Eq[0] * x1)) / Eq[1];
                    y2 = (((-1) * Eq[2]) - (Eq[0] * x2)) / Eq[1];


                    Graphics lines = pictureBox1.CreateGraphics();
                    lines.DrawLine(Pens.Yellow, x1, Convert.ToInt32(y1), x2, Convert.ToInt32(y2));
                
                i = i + 1;
            }
        }

        private void button4_Click(object sender, EventArgs e)
        {
            i = 0;
            listBox1.Items.Clear();

            gImage1 = (Image)defaultImage1.Clone();

            Graphics g = Graphics.FromImage(gImage1);
            pictureBox1.Image = gImage1;

            gImage2 = (Image)defaultImage2.Clone();

            g = Graphics.FromImage(gImage2);
            pictureBox2.Image = gImage2;
        }

        private void button1_Click(object sender, EventArgs e)
        {
            double[,] A = new double[i/2,9];
            double[] w;
            double[,] v;

            int j;
            for ( j =0; j<i/2; j++)
            {
                A[j, 0] = P[j, 0] * PD[j, 0];
                A[j, 1] = P[j, 1] * PD[j, 0];
                A[j, 2] = PD[j, 0];
                A[j, 3] = P[j, 0] * PD[j, 1];
                A[j, 4] = P[j, 1] * PD[j, 1];
                A[j, 5] = PD[j, 1];
                A[j, 6] = P[j, 0];
                A[j, 7] = P[j, 1];
                A[j, 8] = 1;

            }


            svdcmp(A, out w, out v);

            foreach (var k in w)
            {
                Console.WriteLine(k.ToString());
            }

            double min = w[0];
            int flag = 0;
            for (j = 1; j < 9; j++)
            {
                if (w[j] < min)
                {
                    min = w[j];
                    flag = j;
                }
            }

            

            F[0, 0] = v[0, flag];
            F[0, 1] = v[1, flag];
            F[0, 2] = v[2, flag];
            F[1, 0] = v[3, flag];
            F[1, 1] = v[4, flag];
            F[1, 2] = v[5, flag];
            F[2, 0] = v[6, flag];
            F[2, 1] = v[7, flag];
            F[2, 2] = v[8, flag];

            int l;
            Console.WriteLine("");
            for (j = 0; j < 3; j++)
            {
                
                for (l = 0; l < 3; l++)
                {
                    Console.Write("\t" + F[j,l]);
                }
                Console.WriteLine("");
            }

            textBox1.Visible = true;
            textBox1.AppendText(Environment.NewLine);
            textBox1.AppendText(" | Fundamental Matrix | ");
            textBox1.AppendText(Environment.NewLine);
            textBox1.AppendText(Environment.NewLine);
            textBox1.AppendText("[\t" + F[0, 0] + " ,\t" + F[0, 1] + " ,\t" + F[0, 2] + "\t]");
            textBox1.AppendText(Environment.NewLine);
            textBox1.AppendText(Environment.NewLine);
            textBox1.AppendText("[\t" + F[1, 0] + " ,\t" + F[1, 1] + " ,\t" + F[1, 2] + "\t]");
            textBox1.AppendText(Environment.NewLine);
            textBox1.AppendText(Environment.NewLine);
            textBox1.AppendText("[\t" + F[2, 0] + " ,\t" + F[2, 1] + " ,\t" + F[2, 2] + "\t]");
            textBox1.AppendText(Environment.NewLine);
            textBox1.AppendText(Environment.NewLine);
        }

        private void button5_Click(object sender, EventArgs e)
        {

            i = 0;
            listBox1.Items.Clear();

            gImage1 = (Image)defaultImage1.Clone();

            Graphics g = Graphics.FromImage(gImage1);
            pictureBox1.Image = gImage1;

            gImage2 = (Image)defaultImage2.Clone();

            g = Graphics.FromImage(gImage2);
            pictureBox2.Image = gImage2;

            if (epipolar)
            {
                button5.BackColor = Color.Red;
                epipolar = false;
            }
            else
            {
                button5.BackColor = Color.Green;
                epipolar = true;
            }
        }
    }   
}
