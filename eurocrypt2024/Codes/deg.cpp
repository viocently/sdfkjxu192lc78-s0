#include<iostream>
#include<vector>
#include<algorithm>
#include<bitset>

using namespace std;

#define A 0
#define B 93
#define C 177

#define RA 92
#define RB 83
#define RC 110

/*
 * la = S242 + S287 + S68
 * lb = S65 + S92 + S170
 * lc = S161 + S176 + S263
*/

vector<int> K(128);
vector<int> V(128);
vector<int *> D;

const int inf = -100000;

int max( int d1, int d2, int d3 )
{
    int m = -10000000;
    if ( m < d1 ) m = d1;
    if ( m < d2 ) m = d2;
    if ( m < d3 ) m = d3;
    return m;
}

int maxx( int d1, int d2, int d3, int d4 )
{
    int m = -10000000;
    if ( m < d1 ) m = d1;
    if ( m < d2 ) m = d2;
    if ( m < d3 ) m = d3;
    if ( m < d4 ) m = d4;
    return m;
}

int max2( int d1, int d2 )
{
    if ( d1 > d2 ) return d1;
    else return d2;
}


int min( int d1, int d2, int d3 )
{
    int m = 10000000;
    if ( m > d1 ) m = d1;
    if ( m > d2 ) m = d2;
    if ( m > d3 ) m = d3;
    return m;
}

int max7( int d1, int d2, int d3, int d4, int d5, int d6, int d7 )
{
    int m = -10000000;
    if ( m < d1 ) m = d1;
    if ( m < d2 ) m = d2;
    if ( m < d3 ) m = d3;
    if ( m < d4 ) m = d4;
    if ( m < d5 ) m = d5;
    if ( m < d6 ) m = d6;
    if ( m < d7 ) m = d7;
    return m;
}

int At( int t1 )
{
    if ( t1 < 0 )
        return D[0][A - t1];
    else
        return D[t1][A];
}

int Bt( int t1 )
{
    if ( t1 < 0 )
        return D[0][B - t1];
    else
        return D[t1][B];
}

int Ct( int t1 )
{
    if ( t1 < 0 )
        return D[0][C - t1];
    else
        return D[t1][C];
}


int * update( int t )
{
    int ga, gb, gc;
    int t1, t2;
    int d1, d2, d3;

    // ga
    t1 = t - RC;

    if ( t1 <= 0 )
        ga = Ct(t1) + Ct(t1+ 1);
    else
    {
        t2 = t1 - RB; 
        //d1 = min( D[t2][B] + D[t1+1][C], D[t2 + 2][B] + D[t2][C], D[t2][B] + D[t2+1][B] + D[t2+2][B] );
        d1 = min( Bt(t2) + Ct(t1+1), Bt(t2+2 ) + Ct(t1), Bt( t2 ) + Bt( t2+1 ) + Bt( t2+2 ) );

        d2 = max( D[t1][161], D[t1][176], D[t1][263] ) + Ct(t1);
            
        d3 = max( D[t1 - 1][161], D[t1 - 1][176], D[t1 - 1][263] ) + Ct(t1+1);

        ga = max( d1, d2, d3 );
    }

    int la = maxx( D[t - 1][242], D[t - 1][287], D[t - 1][68], K[0] );
	// K circle
	int tmppp = K[0];
	for ( int i = 0; i< 127; i++ )
		K[i] = K[i + 1];
	K[127] = tmppp;

    // gb
    t1 = t - RA;
    if ( t1 <= 0 )
        gb = At(t1) + At(t1+ 1);
    else
    {
/*
 * la = S242 + S287 + S68
 * lb = S65 + S92 + S170
 * lc = S161 + S176 + S263
*/
        t2 = t1 - RC; 
        //d1 = min( D[t2][C] + D[t1+1][A], D[t2 + 2][C] + D[t2][A], D[t2][C] + D[t2+1][C] + D[t2+2][C] );
        d1 = min( Ct(t2) + At(t1+1), Ct(t2 + 2) + At(t1), Ct(t2) + Ct(t2+1) + Ct(t2+2) );
        d2 = max( D[t1][242], D[t1][287], D[t1][68] ) + At(t1);
        d3 = max( D[t1 - 1][242], D[t1 - 1][287], D[t1 - 1][68] ) + At(t1+1);
        gb = max( d1, d2, d3 );
    }

    int lb = maxx( D[t - 1][65], D[t - 1][92], D[t - 1][170], V[0] );
	// V circle
	int tmpp = V[0];
	for ( int i = 0; i< 127; i++ )
		V[i] = V[i + 1];
	V[127] = tmpp;
        
    //gc
    t1 = t - RB;
    if ( t1 <= 0 )
        gc = Bt(t1) + Bt(t1+ 1);
    else
    {
        t2 = t1 - RA; 
        //d1 = min( D[t2][A] + D[t1+1][B], D[t2 + 2][A] + D[t2][B], D[t2][A] + D[t2+1][A] + D[t2+2][A] );
        d1 = min( At(t2) + Bt(t1+1), At(t2 + 2) + Bt(t1), At(t2) + At(t2+1) + At(t2+2) );
        d2 = max( D[t1][65], D[t1][92], D[t1][170] ) + Bt(t1);
        d3 = max( D[t1 - 1][65], D[t1 - 1][92], D[t1 - 1][170] ) + Bt(t1+1);
        gc = max( d1, d2, d3 );
    }

    int lc = max( D[t - 1][161], D[t - 1][176], D[t - 1][263] );
/*
 * la = S242 + S287 + S68
 * lb = S65 + S92 + S170
 * lc = S161 + S176 + S263
*/
    int * p = new int[288];

    for ( int i = 0; i < 288; i++ ) 
        if ( ( i != A ) && ( i != B ) && ( i != C ) ) 
            p[i] = D[t-1][ (i + 287) % 288];
    
    p[A] = max2( ga, la );
    p[B] = max2( gb, lb );
    p[C] = max2( gc, lc );

    //for ( int i = 0; i < 288; i++ )
    //  cout << p[i] << '\t';
    //cout << endl;
    return p;
}

int computeDegree( const bitset<128>& cube, const int rounds, const bitset<544> & intervector  )
//int computeDegree( const bitset<128>& cube, const int rounds )
{
    D.clear();
    int * p = new int[288];

	// key
    for ( int i = 0; i < 93; i++ )
		p[i] = 0;

	for ( int i = 93; i < 93 + 128; i++ )
		if ( cube[i - 93] == 1 )
			p[i] = 1;
		else
			p[i] = inf;

	for ( int i = 93 + 128; i < 287; i++ )
		p[i] = 0;

	p[287] = inf;

    D.push_back( p );

	// initialize K and V
	for ( int i = 0; i < 128; i++ )
		K[i] = 0;

	for ( int i = 0; i < 128; i++ )
		if ( cube[i] == 1 )
			V[i] = 1;
		else
			V[i] = inf;

    int degree = 0;
    for ( int i = 1; i < rounds + 1; i++ )
    {
        D.push_back ( update( i ) );
    }

    for ( int j = 0; j < 128; j++ )
        if ( intervector[j + 288 + 128] == 1 )
            degree += K[j];
    for ( int j = 0; j < 128; j++ )
        if ( intervector[j + 288] == 1 )
            degree += V[j];
    for ( int j = 0; j < 288; j++ )
        if ( intervector[j] == 1 )
            degree += D[rounds][j];

    //auto degree = max7 ( D[rounds][65], D[rounds][92], D[rounds][161],
    //			D[rounds][176], D[rounds][242], D[rounds][287], K[0] ) ; 
    //cout << "ROUND: [" << i << "] " << "DEGREE: [" << deg << "]" << endl;

	for ( auto & it : D )
		delete [] it;

    return degree;
}


/*
int main()
{
    //int I_liu[] = {0, 2, 4, 6, 8, 10, 12, 15, 17, 19, 21, 23, 25, 27, 30, 32, 34, 36, 38, 40, 42, 45, 47, 49, 51, 53, 55, 57, 60, 62, 64, 66, 68, 70, 72, 75, 79};
    vector<int> cube {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 29, 31,
		33, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 57, 59, 61, 63, 65, 67, 69,
		72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 102, 104,
		107, 109, 111, 113, 115, 117, 119, 122, 124, 126 };

	bitset<128> I;
	for ( int i = 0; i < 128; i++ )
		I[i] = 1;

	cout << cube.size() << endl;
	for ( int i = 0; i < 1000; i++ )
		cout << i << " " << computeDegree( I, i ) << endl;
}
*/

