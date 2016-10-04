#include "fft.h"
#include "mtxutl.h"

#define SWAP( a, b ) tempr=( a ); ( a ) = ( b ); ( b ) = tempr;

void four1( float *data, int nn, int isign )
{
	unsigned long n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta;
	float tempr, tempi;

	n = nn << 1;
	j = 1;
	for( i=1; i<n; i+=2 )
	{
		if( j > i )
		{
			SWAP( data[j], data[i] );
			SWAP( data[j+1], data[i+1] );
		}
		m = n >> 1;
		while( m >= 2 && j > m )
		{
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax = 2;
	while( n > mmax )
	{
		istep = mmax << 1;
		theta = isign * ( 2 * PI / mmax );
		wtemp = sin( 0.5 * theta );
		wpr = -2.0 * wtemp * wtemp;
		wpi = sin( theta );
		wr = 1.0;
		wi = 0.0;
		for( m=1; m<mmax; m+=2 )
		{
			for( i=m; i<=n; i+=istep )
			{
				j = i + mmax;
				tempr = wr * data[j] - wi * data[j+1];
				tempi = wr * data[j+1] + wi * data[j];
				data[j] = data[i] - tempr;
				data[j+1] = data[i+1] - tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr = ( wtemp = wr ) * wpr - wi * wpi + wr;
			wi = wi * wpr + wtemp * wpi + wi;
		}
		mmax = istep;
	}
}

int fft( int n, Fukusosuu *in, int disp )
{
	int i, m, sign;
	float *x;
	m = abs( n );
	sign = n / m;
	x = calloc( m*2+1, sizeof( float ) );

	for( i=0; i<m; i++ )
	{
		x[i*2+1] = in[i].R;
		x[i*2+2] = in[i].I;
	}
#if 0
if( disp )                                                      
{                                                                    
    FILE *fp;        
    fp = fopen( "inputOfFft", "w" );                                           
    for( i=0; i<m; i++ )                
        fprintf( fp, "%f %f\n", x[i*2+1], x[i*2+2] );        
    fclose( fp );        
    system( "vi inputOfFft < /dev/tty > /dev/tty " );        
}                
#endif        
	four1( x, m, sign );
#if 0
if( disp )                                                      
{                                                                    
    FILE *fp;        
    fp = fopen( "outputOfFft", "w" );                                           
    for( i=0; i<m; i++ )                
        fprintf( fp, "%f %f\n", x[i*2+1], x[i*2+2] );        
    fclose( fp );        
    system( "vi outputOfFft < /dev/tty > /dev/tty " );        
}                
#endif        
	for( i=0; i<m; i++ )
	{
		in[i].R = x[i*2+1];
		in[i].I = x[i*2+2];
	}
}
		

#if 0


int fft( int n, Fukusosuu *in, int disp )
{
	int i, m;
	float *x, *y;

	m = abs( n );

	x = calloc( 100000, sizeof( float ) );
	y = calloc( 100000, sizeof( float ) );

	for( i=0; i<m; i++ )
	{
		x[i] = (float)in[i].R;
		y[i] = (float)in[i].I;
	}	

	fft_hontai( n, x, y );
	for( i=0; i<m; i++ )
	{
		(in+i)->R = (double)x[i];
		(in+i)->I = (double)y[i];
	}	
	free( x );
	free( y );
}
#endif
