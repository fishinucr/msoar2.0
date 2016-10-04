#include "fft.h"
#include "mtxutl.h"

FILE *gfp;
/*
  {\tt fft()}.
*/
static void make_sintbl(int n, double *sintbl)
{
        int i, n2, n4, n8;
        double c, s, dc, ds, t;

        n2 = n / 2;  n4 = n / 4;  n8 = n / 8;
#if 0
        t = sin(PI / n);
        dc = 2 * t * t;  ds = sqrt(dc * (2 - dc));
        t = 2 * dc; c = 1.0; s = 0.0; 
        sintbl[n4] = 1.0;  
        sintbl[0] = 0.0;
        for (i = 1; i < n8; i++) {
                c -= dc;  dc += t * c;
                s += ds;  ds -= t * s;
                sintbl[i] = s;  sintbl[n4 - i] = c;
        }
        if (n8 != 0) sintbl[n8] = sqrt(0.5);
#else
        sintbl[n4] = 1.0;  
        sintbl[0] = 0.0;
        for (i = 1; i < n4; i++)
			sintbl[i] = sin( (double)i / n * ( 2 * PI ) );
#endif
        for (i = 0; i < n4; i++)
                sintbl[n2 - i] = sintbl[i];
        for (i = 0; i < n2 + n4; i++)
                sintbl[i + n2] = - sintbl[i];
#if 0
{
    FILE *fp;
    fp = fopen( "sintbl", "w" );
    for( i=0; i<n; i++ )
        fprintf( fp, "%d %f\n", i, sintbl[i] );
    fclose( fp );
    fp = fopen( "plot", "w" );
    fprintf( fp, "plot 'sintbl'\npause 5" );
    fclose( fp );
    system( "gnuplot plot" );
}
#endif
}
/*
  {\tt fft()}.
*/
static void make_bitrev(int n, int *bitrev)
{
        int i, j, k, n2;

        n2 = n / 2;  i = j = 0;
        for ( ; ; ) {
                bitrev[i] = j;
                if (++i >= n) break;
                k = n2;
                while (k <= j) {  j -= k;  k /= 2;  }
                j += k;
        }
}
/*
*/
#if 0
int fft(int n, Fukusosuu *x)
{
        static int    last_n = 0;    /*  {\tt n} */
        static int   *bitrev = NULL; /*  */
        static float *sintbl = NULL; /*  */
        int i, j, k, ik, h, d, k2, n4, inverse;
        float t, s, c, dR, dI;

        /*  */
        if (n < 0) {
                n = -n;  inverse = 1;  /*  */
        } else inverse = 0;
        n4 = n / 4;
        if (n != last_n || n == 0) {
                last_n = n;
                if (sintbl != NULL) free((void *)sintbl);
                if (bitrev != NULL) free((void *)bitrev);
                if (n == 0) return 0;  /*  */
                sintbl = malloc((n + n4) * sizeof(float));
                bitrev = malloc(n * sizeof(int));
                if (sintbl == NULL || bitrev == NULL) {
                        fprintf(stderr, "\n");  return 1;
                }
                make_sintbl(n, sintbl);
                make_bitrev(n, bitrev);
        }
        for (i = 0; i < n; i++) {    /*  */
                j = bitrev[i];
                if (i < j) {
                        t = x[i].R;  x[i].R = x[j].R;  x[j].R = t;
                        t = x[i].I;  x[i].I = x[j].I;  x[j].I = t;
                }
        }
        for (k = 1; k < n; k = k2) {    /*  */
                h = 0;  k2 = k + k;  d = n / k2;
                for (j = 0; j < k; j++) {
                        c = sintbl[h + n4];
                        if (inverse) s = - sintbl[h];
                        else         s =   sintbl[h];
                        for (i = j; i < n; i += k2) {
                                ik = i + k;
                                dR = s * x[ik].I + c * x[ik].R;
                                dI = c * x[ik].I - s * x[ik].R;
                                x[ik].R = x[i].R - dR;  x[i].R += dR;
                                x[ik].I = x[i].I - dI;  x[i].I += dI;
                        }
                        h += d;
                }
        }
        if (! inverse)    /* n */
                for (i = 0; i < n; i++) {  x[i].R /= (double)n;  x[i].I /= (double)n;  }
        return 0;  /*  */
}
#else
int fft(int n, Fukusosuu *x, int disp)
{
        static int    last_n = 0;    /*  {\tt n} */
        static int    *bitrev = NULL; /*  */
        static double *sintbl = NULL; /*  */
        int i, j, k, ik, h, d, k2, n4, inverse;
        double t, s, c, dR, dI;
		double tmp1, tmp2;

        /*  */
        if (n < 0) {
                n = -n;  inverse = 1;  /*  */
        } else inverse = 0;
#if 0
if( disp )
{
    FILE *fp;
    fp = fopen( "inputOfFft", "w" );
    for( i=0; i<n; i++ )
        fprintf( fp, "%f %f\n", x[i].R, x[i].I );
    fclose( fp );
    system( "vi inputOfFft < /dev/tty > /dev/tty " );
}
#endif
        n4 = n / 4;
        if (n != last_n || n == 0) {
                last_n = n;
                if (sintbl != NULL) free((void *)sintbl);
                if (bitrev != NULL) free((void *)bitrev);
                if (n == 0) return 0;  /*  */
                sintbl = (double *)malloc((n + n4) * sizeof(double));
                bitrev = (int    *)malloc(n * sizeof(int));
                if (sintbl == NULL || bitrev == NULL) {
                        fprintf(stderr, "\n");  return 1;
                }
                make_sintbl(n, sintbl);
                make_bitrev(n, bitrev);
        }
        for (i = 0; i < n; i++) {    /*  */
                j = bitrev[i];
                if (i < j) {
                        t = x[i].R;  x[i].R = x[j].R;  x[j].R = t;
                        t = x[i].I;  x[i].I = x[j].I;  x[j].I = t;
                }
        }
#if 0
if( disp ) gfp = fopen( "x", "w" );
if( disp ) for( i=0; i<n; i++ )
	fprintf( gfp, "i=%d x=%f+i%f\n", i, x[i].R, x[i].I );
#endif
        for (k = 1; k < n; k = k2) {    /*  */
                h = 0;  k2 = k + k;  d = n / k2;
                for (j = 0; j < k; j++) {
                        c = sintbl[h + n4];
                        if (inverse) s = - sintbl[h];
                        else         s =   sintbl[h];
                        for (i = j; i < n; i += k2) {
                                ik = i + k;
								tmp1 = c * x[ik].I; tmp2 = s * x[ik].R;
								dI = tmp1 - tmp2;
								tmp1 = s * x[ik].I; tmp2 = c * x[ik].R;
								dR = tmp1 + tmp2;
#if 0
if( disp ) fprintf( gfp, "B k,j,i,ik=%d,%d,%d,s=%f,c=%f,x[%d]=%f+-i%f, dR=%f ? %f+%f \n", k, j, i, s, c, i, x[i].R, x[i].I, dR, s * x[ik].I, c * x[ik].R );
if( disp ) fprintf( gfp, "B k,j,i,ik=%d,%d,%d,s=%f,c=%f,x[%d]=%f+-i%f, dR=%f\n", k, j, i, s, c, ik, x[ik].R, x[ik].I, dR );
#endif
                                x[ik].R = x[i].R - dR;  x[i].R += dR;
                                x[ik].I = x[i].I - dI;  x[i].I += dI;
#if 0
if( disp ) fprintf( gfp, "A k,j,i=%d,%d,%d,s=%f,c=%f,x[%d]=%f+-i%f, dR=%f\n", k, j, i, s, c, ik, x[ik].R, x[ik].I, dR );
#endif
                        }
                        h += d;
                }
#if 0
if( disp ) fprintf( gfp, "ik=%d x=%f+i%f\n", ik, x[ik].R, x[ik].I );
#endif
        }
        if (! inverse)    /* n */
                for (i = 0; i < n; i++) {  x[i].R /= (double)n;  x[i].I /= (double)n;  }
#if 0
if( disp ) fclose( gfp );
if( disp ) system( "vi x < /dev/tty > /dev/tty " );
#endif
#if 0
if( disp )
{
	FILE *fp;
    fp = fopen( "outputOfFft", "w" );
    for( i=0; i<n; i++ )
        fprintf( fp, "%f %f\n", x[i].R, x[i].I );
    fclose( fp );
    system( "vi outputOfFft < /dev/tty > /dev/tty " );
}
#endif
        return 0;  /*  */
}
#endif
