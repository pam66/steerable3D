package eu.marbilab.imagej;

/**
 * 
 * <h1>Steerable3D</h1>
 * <p>This plugin applies a steerable filter based on 3D Gaussian derivatives
 * following the prescriptions in<br/>
 * Schneider M, Hirsch S, Weber B, Székely G, Menze BH. Med. Image Anal. 2015; 19: 220–249<br/>
 * Filter response is shown as a new stack.></p>
      
 * @author P. Miocchi (MARBILab - Fondazione Santa Lucia)
 * @version 0.7.2
 * 
 * [------before version control-------------
 * v.0.1:
 * - 31-Aug-17: bug corrected in setConvolution_wr  
 * - 25-Sep-17: bug corrected in setConvolution.myThread.run()  
 * v.0.2:
 * - 10-Oct-17: multithreading improved in setConvolution.myThread.run()  
 * v.0.3:
 * - 27-Oct-17: status bar slightly improved during optimal angle search  
 * v.0.4:
 * - 17-Nov-17: BUG fixed in setConvolution  
 * v.0.5:
 * - 18-Nov-17: angle limits given by user  
 * v.0.6:
 * - 22-Dic-17: bug corrected in optimal angle search  
 * v.0.7:
 * - 27-Dic-17: infos on # threads shown in showProgress  
 * - 7-Feb-18: name changed into Steerable3D_ (the serial version is now named Steerable3Dser_) 
 * ----------------------------------------]
 * put on github as 0.7.0
 **/
import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.filter.*;
import java.util.concurrent.atomic.AtomicInteger;

public class Steerable3D_ implements PlugInFilter {
        
  protected static final String VERSION = "0.7.2";
  protected static final String FILTER_ID = "S3D";
  ImagePlus imp;
  protected ImageStack stack;
  protected float[][][][] conv;
  protected float[][] convrot2;
  protected float[] kernel;
  final double SQRTPI2 = 2.50662827463;
  int[][] bincoeff;
  //---- user parameters
  int M, A, B; // 3d filter order with M >= A >= B 
  int SIGMA; // filter length scale
  double THETA, PSI; //filter orientation
  double thstep, psistep; //filter orientation steps
  double thfrom, thto, psifrom, psito;
  //--------------------------------

  @Override
  public int setup(String arg, ImagePlus imp) {

    if (arg.equals("about")) {
      IJ.showMessage("About Steerable3D...",
              "This plugin applies a steerable filter based on 3D Gaussian derivatives.\n"
              + "Filter feature is shown in a new stack.\n"
              + "Ver. "+VERSION+" by P. Miocchi (MARBILab - Fondazione Santa Lucia).");
      return DONE;
    }
    int dummy, ii, choice;
    boolean template;

    this.bincoeff = new int[5][5];
    this.imp = imp;
    GenericDialog gd = new GenericDialog("Steerable3Dp");
    gd.enableYesNoCancel("Run", "Draw template");
    gd.addHelp("https://imagej.nih.gov/ij/plugins/index.html");
    String[] azione = {"Custom (M,a,b)", "Edge (1,0,0)", "Ridge (2,1,0)", "Vessel (2,0,0)"};
//    ImagePlus logo = new ImagePlus("Steerable3Dp_.tif");    
//    gd.addImage(logo);
    gd.addChoice("Filter type: ", azione, "Custom (M,a,b)");
    //gd.addCheckbox("Draw Template only",false);
    gd.addNumericField("M (<5): ", 2, 0); //M
    gd.addNumericField("a (<= M): ", 0, 0); //A
    gd.addNumericField("b (<= a): ", 0, 0); //B
    gd.addNumericField("scale (pixels): ", 1, 0); //SIGMA
    gd.addNumericField("Elevation step (deg): ", 15., 0); //theta step
    gd.addNumericField("Azimuth step (deg): ", 30., 0); //psi step
    gd.addNumericField("Elevation from (deg): ", 0., 0); //theta from
    gd.addNumericField("Elevation to (deg): ", 90., 0); //theta to
    gd.addNumericField("Azimuth from (deg): ", 0., 0); //psi from
    gd.addNumericField("Azimuth to (deg): ", 360., 0); //psi to
    gd.showDialog();
    if (gd.wasCanceled()) {
      return DONE;
    }
    choice = (int) gd.getNextChoiceIndex();
    template = (boolean) !gd.wasOKed();
    switch (choice) {
      case 0:
        this.M = (int) gd.getNextNumber();
        this.A = (int) gd.getNextNumber();
        this.B = (int) gd.getNextNumber();
        break;
      case 1: //edge
        this.M = 1;
        this.A = 0;
        this.B = 0;
        dummy = (int) gd.getNextNumber();
        dummy = (int) gd.getNextNumber();
        dummy = (int) gd.getNextNumber();
        break;
      case 2: //ridge
        this.M = 2;
        this.A = 1;
        this.B = 0;
        dummy = (int) gd.getNextNumber();
        dummy = (int) gd.getNextNumber();
        dummy = (int) gd.getNextNumber();
        break;
      case 3: //vessel
        this.M = 2;
        this.A = 0;
        this.B = 0;
        dummy = (int) gd.getNextNumber();
        dummy = (int) gd.getNextNumber();
        dummy = (int) gd.getNextNumber();
        break;
    }
    this.SIGMA = (int) gd.getNextNumber();
    if (template) {
      this.SIGMA = 5;
    } else {
      try {
        this.stack = imp.getStack();
      } catch (NullPointerException e) {
        IJ.error("Stack required!");
        return DONE;
      }
      this.conv = new float[this.M + 1][this.M + 1][this.stack.getSize()][this.stack.getWidth() * this.stack.getHeight()];
      this.convrot2 = new float[this.stack.getSize()][this.stack.getWidth() * this.stack.getHeight()];
    }
    if (this.M > 4 || this.M < 0 || this.B > this.A || this.A > this.M || this.SIGMA < 1) {
      IJ.error("Wrong parameters!\nCheck limits and constraints.");
      return DONE;
    }
    this.thstep = (double) gd.getNextNumber();
    this.psistep = (double) gd.getNextNumber();

    this.thfrom = (double) gd.getNextNumber();
    this.thto = (double) gd.getNextNumber();

    this.psifrom = (double) gd.getNextNumber();
    this.psito = (double) gd.getNextNumber();

    if (this.psifrom > this.psito || this.psifrom > 360. || this.psifrom < 0.
            || this.psito > 360. || this.psito < 0.
            || this.thfrom > this.thto || this.thfrom > 90. || this.thfrom < 0.
            || this.thto > 90. || this.thto < 0.) {
      IJ.error("Wrong angle values!");
      return DONE;
    }
    //this.M=1; this.A=0; this.B=0; this.SIGMA=2;
    this.THETA = 0.;
    this.PSI = 0.; //default: along x on the xy plane
    final int kernelRadius = 3 * this.SIGMA;
    this.kernel = new float[(2 * kernelRadius + 1)*(2 * kernelRadius + 1)*(2 * kernelRadius + 1)];

    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 5; j++) {
        this.bincoeff[i][j] = this.BinCoeff(i, j);
      }
    }

    if (template) {
      //--- create an image of a kernel
      ImagePlus kimage = NewImage.createFloatImage("kernel ("
              + Integer.toString(this.M) + "," + Integer.toString(this.A) + "," + Integer.toString(this.B)
              + ")", 2 * kernelRadius + 1, 2 * kernelRadius + 1, 2 * kernelRadius + 1, NewImage.GRAY32);
      ImageProcessor newip = kimage.getProcessor();
      ImageStack kstack = kimage.getStack();
      float[] kpixels;
      this.setKernel(this.M, this.A, this.B, kernelRadius);
      final int stride=2 * kernelRadius + 1;
      for (int k = 1; k <= stride; k++) {
        ii = 0;
        kpixels = (float[]) kstack.getPixels(k);
        for (int j = 0; j < stride; j++) {
          for (int i = 0; i < stride; i++) {
            kpixels[ii] = this.kernel[(k - 1)*stride*stride+j*stride+i];
            ii++;
          }
        }
      }
      newip.resetMinAndMax();
      kimage.show();
      kimage.getCanvas().setMagnification(20. / (double) this.SIGMA);
      kimage.updateAndDraw();
      return DONE;
    }
    //------------

    return DOES_32 + STACK_REQUIRED + NO_CHANGES;
  }

  @Override
  public void run(ImageProcessor ip) {
    int ic, jc, kc, ii;
    final int kernelRadius = 3 * this.SIGMA;
    ImagePlus newimage = this.imp.duplicate();
    ImageStack newstack = newimage.getStack();
    Thread[] threads;
    //this is just to find the max numbers of PEs    
    threads = ij.util.ThreadUtil.createThreadArray();
    final int npmax = threads.length;
    int np = this.stack.getHeight() * this.stack.getWidth() * this.stack.getSize();
    //determine np as the optimal no. of PEs such that each PE has <= 1e7 voxels 
    np = (int) (np / 1e7); 
    if (np < 1) {
      np = 1;
    }
    if (np > npmax) {
      np = npmax;
    }
    //reallocate np threads
    threads = ij.util.ThreadUtil.createThreadArray(np);

    int ntot = (this.M + 1) * (this.M + 2) / 2;
    int n = 0;
    for (int i = 0; i <= this.M; i++) {
      for (int j = 0; j <= i; j++) {
        IJ.showStatus("Calculating convolution w/o rotation..." + " (#threads = " + Integer.toString(npmax) + ")");
        this.setConvolution_wr(this.M, i, j, kernelRadius, n, ntot);
        n += this.stack.getSize() * this.stack.getHeight();
      }
    }
    IJ.showStatus("done.");

    for (int k = 0; k < this.stack.getSize(); k++) {
      for (ii = 0; ii < this.stack.getWidth() * this.stack.getHeight(); ii++) {
        this.convrot2[k][ii] = -Float.MAX_VALUE;
      }
    }

    ntot = ((int) ((this.thto - this.thfrom) / this.thstep) + 1)
            * ((int) ((this.psito - this.psifrom) / this.psistep) + 1);
    n = 0;
    IJ.showProgress(n, ntot);
    for (this.THETA = this.thfrom; this.THETA <= this.thto; this.THETA += this.thstep) {
      for (this.PSI = this.psifrom; this.PSI <= this.psito; this.PSI += this.psistep) {
        if (this.PSI == 360.) {
          continue;
        }
        IJ.showStatus("Applying rotation to find the optimal angle..." + " (#threads = " + Integer.toString(threads.length) + ")");
        this.setConvolution(this.M, this.A, this.B, Math.toRadians(this.THETA), Math.toRadians(this.PSI), threads);
        n++;
        IJ.showProgress(n, ntot);
      }
    }

    IJ.showStatus("done.");

    for (int k = 0; k < this.stack.getSize(); k++) {
      newstack.setPixels(this.convrot2[k], k + 1);
    }
    newimage.setTitle(imp.getTitle()+"_"+FILTER_ID+ "_Response("
            + Integer.toString(this.M) + "," + Integer.toString(this.A) + "," + Integer.toString(this.B) +
            "," + Integer.toString(this.SIGMA)
            + ")");
    newimage.updateAndDraw();
    ImageProcessor newip2 = newimage.getProcessor();
    newip2.resetMinAndMax();
    newimage.show();

  }

  /*
   evaluate the kernel, given {m,a,b}
   */
  protected void setKernel(int m, int a, int b, int kernelRadius) {
    final int stride=2 * kernelRadius + 1;
    int addr;
    for (int kc = 0; kc < stride; kc++) {
      addr=kc*stride*stride;
      for (int jc = 0; jc < stride; jc++) {
        for (int ic = 0; ic < stride; ic++) {
          this.kernel[addr+ic] = (float) this.dGauss3D(ic - kernelRadius, jc - kernelRadius, kc - kernelRadius, this.SIGMA, m, a, b);
        }
        addr=addr+stride;
      }
    }
  }

  /*
   do convolution WITH rotation of a STF of order m,a,b (eq. 6, Schneider et al, 2015, MIA, 220)
   and put it into <convrot2>.
   setConvolution_wr must be called first to evaluate convolution without rotation
   */
  protected void setConvolution(final int m, final int a, final int b, final double theta, final double psi,
          Thread[] threads) {

    final ImageStack sstack = this.stack;
    final float[][] convrot = new float[this.stack.getSize()][this.stack.getWidth() * this.stack.getHeight()];
    final float[][] sconvrot2 = this.convrot2;
    final float[][][][] sconv = this.conv;
    final double sinth = Math.sin(theta);
    final double costh = Math.cos(theta);
    final double sinps = Math.sin(psi);
    final double cosps = Math.cos(psi);
    final int[][] bbincoeff = this.bincoeff;
    class myThread extends Thread {

      private int myid;
      private int threadsNo;

      public myThread(final int _myid, final int _threadsNo) {
        setPriority(Thread.NORM_PRIORITY);
        this.myid = _myid;
        this.threadsNo = _threadsNo;
      }

      @Override
      public void run() {
        // if (k == 0) this.first = true;
        final int quant = sstack.getSize() / this.threadsNo;
        final int from = quant * this.myid;
        final int to = this.myid == this.threadsNo - 1 ? sstack.getSize() : quant * (this.myid + 1);
        float w;
        double wt;
        int i, j, k, ii;

        for (i = 0; i <= m; i++) {
          for (j = 0; j <= i; j++) {


            //----evaluate w^{ij}_{mab}
            w = (float) 0.;

            for (int u1 = 0; u1 <= m - a; u1++) {
              for (int v1 = 0; v1 <= a - b; v1++) {
                for (int w1 = 0; w1 <= b; w1++) {
                  if (u1 + v1 + w1 != i) {
                    continue;
                  }
                  for (int u2 = 0; u2 <= u1; u2++) {
                    for (int w2 = 0; w2 <= w1; w2++) {
                      if (u2 + w2 != j) {
                        continue;
                      }

                      wt =
                              (double) (bbincoeff[m - a][u1] * bbincoeff[a - b][v1] * bbincoeff[b][w1] * bbincoeff[u1][u2]
                              * bbincoeff[w1][w2])
                              * Math.pow(costh, m - a - u2 + w2) * Math.pow(cosps, m - a + b - u1 + v1 - w1)
                              * Math.pow(sinth, b + u2 - w2) * Math.pow(sinps, a - b + u1 - v1 + w1 - u2 - w2);
                      if ((a - v1 - w2) % 2 == 1) {
                        wt = -wt;
                      }
                      w += (float) wt;

                    }
                  }
                }
              }
            }
            for (k = from; k < to; k++) {
              for (ii = 0; ii < sstack.getWidth() * sstack.getHeight(); ii++) {
                convrot[k][ii] += w * sconv[i][j][k][ii];
              }
            }

          }
        }
        for (k = from; k < to; k++) {
          for (ii = 0; ii < sstack.getWidth() * sstack.getHeight(); ii++) {
            sconvrot2[k][ii] = Math.max(sconvrot2[k][ii], convrot[k][ii]);
          }
        }

      }
    }


    for (int ithread = 0; ithread < threads.length; ithread++) {
      threads[ithread] = new myThread(ithread, threads.length);
    }
    ij.util.ThreadUtil.startAndJoin(threads);

  }
  /*
   do convolution without rotation of a STF of order m,a,b with the given
   kernel radius.
   At the end <conv> is equal to f_{mab}^\sigma(I,x) (eq.6, Schneider et al, 2015, MIA, 220)
   */

  protected void setConvolution_wr(int m, final int a, final int b, final int kernelRadius, final int n, final int ntot) {


    // set the kernel depending on m,a,b
    this.setKernel(m, a, b, kernelRadius);
    // do convolution with the 3D image (the stack)

    final Thread[] threads = ij.util.ThreadUtil.createThreadArray();
    final AtomicInteger ai = new AtomicInteger(0);
    final ImageStack sstack = this.stack;
    final float[][][][] sconv = this.conv;
    final float[] skernel = this.kernel;

    class myThread extends Thread {

      private boolean master;

      public myThread(int _myid) {
        setPriority(Thread.NORM_PRIORITY);
        this.master = _myid == 0;
      }

      @Override
      public void run() {
        float[] source;
        int ii;
        final int stride=2 * kernelRadius + 1;

        for (int k = ai.getAndIncrement(); k < sstack.getSize(); k = ai.getAndIncrement()) {
          ii = 0; //k+1-th slice
          for (int j = 0; j < sstack.getHeight(); j++) {
            for (int i = 0; i < sstack.getWidth(); i++) {
              //it is evaluating the (i,j,k) voxel
              sconv[a][b][k][ii] = (float) 0.;
              //convolve the 3D image with the kernel around (i,j,k)
              for (int kc = k - kernelRadius; kc <= k + kernelRadius; kc++) {
                if (kc < 0 || kc >= sstack.getSize()) {
                  continue;
                }
                source = (float[]) sstack.getPixels(kc + 1);
                for (int jc = j - kernelRadius; jc <= j + kernelRadius; jc++) {
                  if (jc < 0 || jc >= sstack.getHeight()) {
                    continue;
                  }
                  for (int ic = i - kernelRadius; ic <= i + kernelRadius; ic++) {
                    if (ic < 0 || ic >= sstack.getWidth()) {
                      continue;
                    }
                    sconv[a][b][k][ii] += source[ic + jc * sstack.getWidth()]
                            * skernel[
                      (kc - k + kernelRadius)*stride*stride+(jc - j + kernelRadius)*stride+ic - i + kernelRadius];
                  }
                }
              }
              ii++;
            }
            if (this.master) {
              IJ.showProgress(n + k * sstack.getHeight() + j, ntot * sstack.getSize() * sstack.getHeight());
            }
          }

        }

      }
    }


    for (int ithread = 0; ithread < threads.length; ithread++) {

      // Concurrently run in as many threads as CPUs  

      threads[ithread] = new myThread(ithread);
    }

    ij.util.ThreadUtil.startAndJoin(threads);

  }

  /*
   gives the partial derivative G^\sigma_{m,a,b}= \delta_x^(m-a) \delta_y^(a-b) \delta_z^b of the 3D
   gaussian with the given sigma. See eq. (3) of Schneider et al, 2015, MIA, 220
   */
  private double dGauss3D(int i, int j, int k, int sigma, int m, int a, int b) {
    return this.dGauss(i, sigma, m - a) * this.dGauss(j, sigma, a - b) * this.dGauss(k, sigma, b);
  }

  /*
   gives the n-th derivative of the gaussian with 0 average and the given sigma
   */
  private double dGauss(int x, int sigma, int n) {
    if (n < 0) {
      return Double.NaN;
    }
    final double s = (double) sigma;
    final double f = 1. / (this.SQRTPI2 * s);
    final double r = ((double) x) / s;

    return (f / Math.exp(r * r * 0.5) * this.Hermite(r, n));
  }

  private double Hermite(double x, int n) {
    switch (n) {
      case 0:
        return 1.;
      case 1:
        return (x);
      case 2:
        return (x * x - 1.0);
      case 3:
        return (Math.pow(x, 3) - 3.0 * x);
      case 4:
        return (Math.pow(x, 4) - 6.0 * x * x + 3.0);
      default:
        return Double.NaN;
    }
  }
  /*
   evaluate a over b, i.e. a!/[b!*(a-b)!]
   */

  private int BinCoeff(int a, int b) {
    if (b > a) {
      return 0;
    }
    if (b == a) {
      return 1;
    }
    int c1 = 1;
    int c2 = 1;
    int i;
    for (i = a - b + 1; i <= a; i++) {
      c1 *= i;
    }
    for (i = 1; i <= b; i++) {
      c2 *= i;
    }
    return (c1 / c2);
  }
}
