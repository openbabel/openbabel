/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.02
 * October 31, 2008
 * Developed at NIST
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC);
 * you can redistribute this software and/or modify it under the terms of 
 * the GNU Lesser General Public License as published by the Free Software 
 * Foundation:
 * http://www.opensource.org/licenses/lgpl-license.php
 */


/* Draw input atom -- Win32 specific */

#include "mode.h"

#ifndef INCHI_ANSI_ONLY /* { */

#ifdef WIN32 /* { */
#include <windows.h>
#include <stdio.h>
#include <math.h>

#include "inpdef.h"
#include "util.h"
#include "dispstru.h"
#include "extr_ct.h"
#include "ichicomp.h"

/* Font size */
#define FONT_NAME "Arial"     /* "MS Sans Serif"; */
/* rgb colors */
#define CLR_BLUE     RGB(  0,   0, 255)
#define CLR_GREEN    RGB(  0, 128,   0)
#define CLR_RED      RGB(255,   0,   0)
#define CLR_PINK     RGB(255, 128, 128)
#define CLR_CYAN     RGB(  0, 255, 255)
#define CLR_LTGREEN  RGB(128 ,255, 128)
#define CLR_LTPURPLE RGB(255 ,0xCC, 255)
#define CLR_YELLOW   RGB(255, 255,   0)
#define CLR_BLACK    RGB(  0,   0,   0)
#define CLR_LTGRAY   RGB(0xCC, 0xCC, 0xCC)
#define CLR_MAGENTA  RGB(255,0,255)
#define CLR_WHITE    RGB(255, 255, 255)

/* local prototypes */
HWND GetConsoleHwnd(void);

#define MY_TIMER_ID 1

typedef struct Box {
    int             xhigh, xlow;
    int             yhigh, ylow;
}                BOX;


/* local prototypes */
int         DrawBond( HDC pDC, int x1, int y1, int x2, int y2, int b_type, int b_stereo, int b_parity, int bInvertBonds, COLORREF clrPen, int nPenWidth );
int         DrawBondStereo( HDC pDC, int x1, int y1, int x2, int y2, int b_stereo, int b_highlight, int bInvertBonds, COLORREF clrPen, int nPenWidth );
int         DrawBondNoStereo( HDC pDC, int x1, int y1, int x2, int y2, int b_type, int b_highlight, COLORREF clrPen, int nPenWidth );
int         DrawBondParity( HDC pDC, int x1, int y1, int x2, int y2, int parity_mark );
void        DrawLine( HDC, int, int, int, int );
void        DrawPenColorFilledPolygon( HDC pDC, const POINT* pnt, int num );
int         DrawTextColorDot( HDC pDC );
int         DrawString( HDC pDC, char *st1, int shift, int x, int y );
int         DrawPreparedString( HDC pDC, char *st1, int shift, int x, int y, int bHighlightTheAtom );
int         DrawColorString( HDC pDC, const char *st, int xs, int ys, int bHighlightTheAtom );
            /* structure drawing */
int         DrawStructure( HDC pDC, inp_ATOM *at, INF_ATOM_DATA *inf_at_data, int num_at, int xoff, int yoff, COLORREF clrPen, int nPenWidth );
            /* all drawing, including text strings and table */
int         DrawTheInputStructure( inp_ATOM *at, INF_ATOM_DATA *inf_at_data, int num_at,
                              HDC pDC, int tx_off, int ty_off, int xoff, int yoff,
                              int width_pix, int height_pix, int bDraw, int bOrigAtom, COLORREF clrPen, int nPenWidth );
            /* calculate sizes, run drawing */
int         CreateInputStructPicture( HDC hDC, MY_WINDOW_DATA *pWinData, RECT *rc, int bPrint, AT_NUMB nNewEquLabel );

void        FreeWinData( MY_WINDOW_DATA* pWinData );
void        InpStructureMarkEquComponents( MY_WINDOW_DATA *pWinData, AT_NUMB nNewEquLabel,
                                   inp_ATOM *at0, inp_ATOM *at1, inf_ATOM *inf_at, int num_at );
int         MyTextOutABC( const char *p, int iFst, int iLst, HDC pDC );

/* window procedure */
LRESULT CALLBACK WndProcDisplayInputStructure(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);

/* main drawing function: create window, save drawing parameters in the window */
int DisplayInputStructure( char *szOutputString, inp_ATOM  *at, INF_ATOM_DATA  *inf_at_data, int num_at, DRAW_PARMS *dp /*, int bTaut, unsigned long ulDisplTime, long *rcPict, int nFontSize*/ );


int         GetFontHeight( HDC );
int         GetFontAscent( HDC );
int         GetFontDescent( HDC pDC );
int         GetFontAveWidth( HDC );
int         GetStringWidth( HDC  pDC, char *pString );
int         GetOneCharInStringWidth( HDC  pDC, const char *pString );
int         MoveHydrogenAtomToTheLeft( char *s, int start, int H );

int         nRound( double X );
double      RoundDouble(double X);
void        roundoff_coord( double dx1, double dx2, int *new_ix1, int *new_ix2 );
BOOL        ReallySetForegroundWindow(HWND hWnd);
HWND        GetConsoleHwnd(void);


/*********************************************************************/
char szWindowClassName[] = "INChI_DrawStructWnd";

/*********************************************************************/

void PrintFileName( const char *fmt, FILE *output_file, const char *szFname )
{
    char szBuf[_MAX_PATH];
    long lBufLen = sizeof(szBuf);
    long lReqBufLen;
    char *pName;
    const char *p;

    lReqBufLen = GetFullPathName( szFname, lBufLen, szBuf, &pName );
    
    if ( lReqBufLen && lReqBufLen < lBufLen ) {
        p = szBuf;
    } else {
        p = szFname;
    }
    fprintf( output_file, fmt, p );
}

/*********************************************************************/
BOOL ReallySetForegroundWindow(HWND hWnd)
{
    BOOL  retVal = FALSE;
    HWND  hForegroundWnd;
    if ( hWnd && IsWindow(hWnd) && (hForegroundWnd = GetForegroundWindow()) ) {
        DWORD dwWindowThreadProcessId = GetWindowThreadProcessId(hForegroundWnd, NULL);
        DWORD dwCurrentThreadId       = GetCurrentThreadId();
        AttachThreadInput(dwWindowThreadProcessId, dwCurrentThreadId, TRUE);
        retVal = SetForegroundWindow(hWnd);
        AttachThreadInput(dwWindowThreadProcessId, dwCurrentThreadId, FALSE);
    }
    return retVal;
}
/*********************************************************************/
HWND GetConsoleHwnd(void)
{
#define MY_BUFSIZE 1024 /* Buffer size for console window titles. */
   HWND hwndFound;         /* This is what is returned to the caller. */
   char pszNewWindowTitle[MY_BUFSIZE]; /* Contains fabricated */
                                       /* WindowTitle. */
   char pszOldWindowTitle[MY_BUFSIZE]; /* Contains original */
                                       /* WindowTitle. */

   /* Fetch current window title. */

   GetConsoleTitle(pszOldWindowTitle, MY_BUFSIZE);

   /* Format a "unique" NewWindowTitle. */

   wsprintf(pszNewWindowTitle,"InChITmpWnd%ul/%ul",
               (unsigned long)GetTickCount(),
               (unsigned long)GetCurrentProcessId());

   /* Change current window title. */
   SetConsoleTitle(pszNewWindowTitle);

   /* Ensure window title has been updated. */

   Sleep(40);

   /* Look for NewWindowTitle. */

   hwndFound=FindWindow(NULL, pszNewWindowTitle);

   /* download INChI packages
   ShellExecute(hwndFound, "open", "http://www.iupac.org/inchi", "", "C:\\", SW_SHOWNORMAL);
   */

   /* Restore original window title. */

   SetConsoleTitle(pszOldWindowTitle);

   return(hwndFound);
}
/****************************************************************************/
int         GetFontHeight( HDC pDC )
{
    TEXTMETRIC    TextMetric;

    GetTextMetrics( pDC, &TextMetric );

    return TextMetric.tmHeight;
}
/****************************************************************************/
int         GetFontAscent( HDC pDC )
{
    TEXTMETRIC    TextMetric;

    GetTextMetrics( pDC, &TextMetric );

    return TextMetric.tmAscent;
}
/****************************************************************************/
int         GetFontDescent( HDC pDC )
{
    TEXTMETRIC    TextMetric;

    GetTextMetrics( pDC, &TextMetric );

    return TextMetric.tmDescent;
}

/****************************************************************************/
int         GetFontAveWidth( HDC  pDC)
{
    TEXTMETRIC    TextMetric;

    GetTextMetrics( pDC, &TextMetric );

    return TextMetric.tmAveCharWidth;
}
/****************************************************************************/
int         GetSubstringWidth( HDC  pDC, int len, char *pString)
{
    SIZE    Size;
    int     widthABC, i;
    ABC     abc;
    GetTextExtentPoint32( pDC, pString, len, &Size );
    for ( i = widthABC = 0; i < len; i ++ ) {
        if ( GetCharABCWidths( pDC,   /* handle to DC */
                               (int)pString[i],  /* first character in range */
                               (int)pString[i],  /* last character in range */
                               &abc ) /* array of character widths */
                              ) {
            widthABC += (int)abc.abcB + abs(abc.abcA) + abs(abc.abcC);
        } else {
            break;
        }
    }
    if ( widthABC > Size.cx ) {
        return (widthABC + 2*Size.cx)/3;  /* hunch */
    }
    return    Size.cx;
}
/****************************************************************************/
void       GetTextSize( HDC pDC, int len, char *pString, int *width, int *height )
{
    SIZE    Size;

    GetTextExtentPoint32( pDC, pString, len, &Size );
    *width = Size.cx;
    *height = Size.cy;
}
/****************************************************************************/
void       GetVertTextSize( HDC pDC, int len, char *pString, int *width, int *height )
{
    SIZE    Size;
    TEXTMETRIC    TextMetric;
    int     i;
    GetTextMetrics( pDC, &TextMetric );
    *height = len * TextMetric.tmHeight;
    *width = 0;
    for ( i = 0; i < len; i ++ ) {
        GetTextExtentPoint32( pDC, &pString[i], 1, &Size );
        *width = inchi_max( *width, (int)Size.cx );
    }
}
/****************************************************************************/
BOOL TextOutVert(
  HDC pDC,           /* handle to DC */
  int nXStart,       /* x-coordinate of starting position */
  int nYStart,       /* y-coordinate of starting position */
  LPCTSTR lpString,  /* character string */
  int cbString,      /* number of characters */
  int cell_width     /* width for center alignment */
  ) {
    TEXTMETRIC    TextMetric;
    int i, dy, ret, char_width;
    GetTextMetrics( pDC, &TextMetric );
    dy = TextMetric.tmHeight;
    for ( i = 0, ret = 1; ret && i < cbString; nYStart += dy, i ++ ) {
        char_width = GetOneCharInStringWidth( pDC, lpString+i );
        ret = TextOut( pDC, nXStart+(cell_width-char_width)/2, nYStart, lpString+i, 1 );
    }
    return ret;
}
/****************************************************************************/
BOOL TextOutHoriz(
  HDC pDC,           /* handle to DC */
  int nXStart,       /* x-coordinate of starting position */
  int nYStart,       /* y-coordinate of starting position */
  LPCTSTR lpString,  /* character string */
  int cbString,      /* number of characters */
  int cell_width
  ) {
    int dX = (cell_width && cbString == 1)? (cell_width - GetOneCharInStringWidth( pDC, lpString ))/2:0;
    return TextOut( pDC, nXStart+dX, nYStart, lpString, cbString );
}
/****************************************************************************/
int         GetStringWidth( HDC  pDC, char *pString)
{
    SIZE    Size;

    GetTextExtentPoint32( pDC, pString, strlen( pString ), &Size );
    return    Size.cx;
}

/****************************************************************************/
int         GetOneCharInStringWidth( HDC  pDC, const char *pString)
{
    SIZE    Size;
    GetTextExtentPoint32( pDC, pString, 1, &Size );

    return    Size.cx;
}


/****************************************************************************/
int DrawStructure( HDC pDC, inp_ATOM *at, INF_ATOM_DATA *inf_at_data, int num_at, int xoff, int yoff, COLORREF clrPen, int nPenWidth )
{
    int             i, next, k;
    int             j, r, shift, b_parity=0, bDraw;
    char            str[64], *atname; /*str[sizeof(st->str[0])+1]; */
    inf_ATOM       *inf_at = inf_at_data? inf_at_data->at : NULL;
    int             bUseInvFlags = (!inf_at_data)? 0 : inf_at_data->pStereoFlags? 1 : 2;
    int             bInvertBonds = 0;
    /*
    int             bInvertBonds = inf_at_data && (inf_at_data->StereoFlags & INF_STEREO_INV);
    */
    /* draw all straight lines (bonds) */

    for ( i = 0; i < num_at; i ++ ) {
        /* draw atom #i. */
        for ( j = 0; j < at[i].valence; j++ ) {
            next = at[i].neighbor[j];
            bDraw = 1;
            /* normally draw bonds if next > i; exception: disconnected terminal hydrogen atoms */
            if ( next < i ) {
                /* check if it is a disconnected terminal atom. Disconnected atoms */
                /* have bonds to the rest of the structure; the bond from the rest */
                /* of the structure to the removed terminal atom has been removed. */
                for ( r = 0; r < at[next].valence; r ++ ) {
                    if ( at[next].neighbor[r] == i ) {
                        bDraw = 0;
                        break;
                    }
                }
            }
            if ( bDraw ) {
                /* at[i].bond_stereo[j] is negative if the pointing wedge */
                /* of a stereo bond is at the at[next] atom */
                if ( inf_at ) {
                    for ( k =0, b_parity = 0; k < MAX_STEREO_BONDS && inf_at[i].cStereoBondParity[k]; k ++ ) {
                        if ( inf_at[i].cStereoBondNumber[k] == j ) {
                            b_parity = inf_at[i].cStereoBondParity[k];
                            if ( inf_at[i].cStereoBondWarning[k] ) {
                                b_parity = -b_parity;
                            }
                            break;
                        }
                    }
                }
                switch( bUseInvFlags ) {
                case 1:
                    if ( at[i].component <= inf_at_data->num_components ) {
                        bInvertBonds = (0 != (inf_at_data->pStereoFlags[at[i].component] & INF_STEREO_INV));
                    } else {
                        bInvertBonds = 0;
                    }
                    break;
                case 2:
                    bInvertBonds = 0 != (inf_at_data->StereoFlags & INF_STEREO_INV);
                    break;
                }
                DrawBond( pDC, nRound(at[i].x+xoff), nRound(at[i].y+yoff),
                               nRound(at[next].x+xoff), nRound(at[next].y+yoff),
                               at[i].bond_type[j], at[i].bond_stereo[j], b_parity, bInvertBonds, clrPen, nPenWidth);
            }
        }
    }
    /* write all strings (heteroatoms + H) */
    if ( inf_at ) {
        for ( i = 0; i < num_at; i++ ) {
            strcpy( str, inf_at[i].at_string );
            DrawPreparedString( pDC, str, -inf_at[i].DrawingLabelLeftShift,    nRound(at[i].x+xoff), nRound(at[i].y+yoff), inf_at[i].cHighlightTheAtom );
        }
    } else {
        /* version which does not use inf_at */
        for ( i = 0; i < num_at; i++ ) {

            /* the direction of the shift: */
            shift = at[i].bDrawingLabelLeftShift;
            /* input structure before normalizing: */
            /* terminal H atoms have not been disconnected; */
            /* isotopic H atoms numbers have not been added to num_H */
            atname = at[i].elname;
            j = 0;
            k = 0;
            if ( at[i].iso_atw_diff && (!atname[0] || isupper(UCINT atname[0])) ) {
                int atw = get_atw_from_elnum( (int)at[i].el_number );
                if ( atw ) {
                    k += sprintf( str+k, "^%d", atw+at[i].iso_atw_diff-(at[i].iso_atw_diff>0));
                } else {
                    k += sprintf( str+k, "^+%d", at[i].iso_atw_diff-(at[i].iso_atw_diff>0));
                }
            }
            /* obsolete section, this never happens for now */
            if ( atname[0] && atname[0] < ' ' && atname[1] ) {
                /* special encoding. The 1st byte contains the 1st delimiter; next delimiters are '/'. */
                /* this allows to draw up to 5 numbers from 1..255 range. */
                int t;
                int ch = ' ' + atname[0]; /* delimiter */
                k += sprintf( str+k, "%u", (unsigned)(unsigned char)atname[1]);
                for ( t = 2; t < sizeof(at->elname); t ++ ) {
                    if ( atname[t] ) {
                        if ( t == 3 )
                            ch = ',';  /* comma after 2nd number to separate tautomer group info */
                        if ( ch != '/' && ch != ',' )
                            k += sprintf( str+k, "(%c)%u", ch, (unsigned)(unsigned char)atname[t]);
                        else
                            k += sprintf( str+k, "%c%u", ch, (unsigned)(unsigned char)atname[t]);

                        ch = '/';
                    }
                }
            } else {
                /* this is the main section to display hydrogen atoms, charges, radicals */
                strncpy( str+k, atname+j, sizeof(at->elname)-j );
                str[sizeof(at->elname)+k-j] = '\0';
                k = strlen( str );
                if ( at[i].num_H ) {
                    strcat( str, "H" );
                    k ++;
                    if ( at[i].num_H > 1 ) {
                        k += sprintf( str+k, "%d", (int)at[i].num_H );
                    }
                }
                for ( j = 0; j < NUM_H_ISOTOPES; j ++ ) {
                    if ( at[i].num_iso_H[j] ) {
                        if ( j == 0 || j !=1 && j != 2 ) {
                            k += sprintf( str+k, "^%dH", j+1 );
                        } else {
                            k += sprintf( str+k, j == 1? "D" : "T" );
                        }
                        if ( at[i].num_iso_H[j] > 1 ) {
                            k += sprintf( str+k, "%d", (int)at[i].num_iso_H[j] );
                        }
                    }
                }
                if ( abs(at[i].charge) > 1 )
                    sprintf( str+k, "%+d", at[i].charge );
                else
                if ( abs(at[i].charge) == 1 )
                    strcat( str, at[i].charge>0? "+" : "-" );
                if ( at[i].radical )
                    strcat( str, at[i].radical==RADICAL_SINGLET? ":" :
                                 at[i].radical==RADICAL_DOUBLET? "." :
                                 at[i].radical==RADICAL_TRIPLET? "^^" : "?");
                k = strlen(str);
                sprintf( str+k, "/%d", i+1 ); /* atom ordering number, 1,2,... */
            }

            DrawString( pDC, str, shift,    nRound(at[i].x+xoff), nRound(at[i].y+yoff) );
        }
    }
    return ( 1 );
}

/****************************************************************************/
int DrawBond( HDC pDC, int x1, int y1, int x2, int y2, int b_type, int b_stereo, int b_parity, int bInvertBonds, COLORREF clrPen, int nPenWidth )
{
    int abs_value_of_stereo  = abs( b_stereo );
    int bond_parity_mark     = b_type & BOND_MARK_PARITY;
    int bond_highlight       = b_type & BOND_MARK_HIGHLIGHT;
    int ret;
    b_type ^= (bond_parity_mark | bond_highlight);

    if ( !b_stereo ||
         /* => no stereo */
         b_type== BOND_TYPE_DOUBLE && abs_value_of_stereo!=STEREO_DBLE_EITHER ||
         /* => ignore unknown double bond stereo */
         b_type >= BOND_TYPE_TRIPLE ||
         /* => ignore bond stereo in case of triple bonds */
         b_type==BOND_TYPE_SINGLE && !(abs_value_of_stereo==STEREO_SNGL_UP ||
                                       abs_value_of_stereo==STEREO_SNGL_EITHER ||
                                       abs_value_of_stereo==STEREO_SNGL_DOWN)
         /* => ignore unknown single bond stereo */ ) {
    /*if ( b_type || !abs_value_of_stereo || (abs_value_of_stereo != 1 && abs_value_of_stereo != 4 && abs_value_of_stereo != 6) ) */
        ret = DrawBondNoStereo( pDC, x1, y1, x2, y2, b_type, 0!=bond_highlight, clrPen, nPenWidth );
    } else {
        ret = DrawBondStereo( pDC, x1, y1, x2, y2, b_stereo, 0!=bond_highlight, bInvertBonds, clrPen, nPenWidth );
    }
    if ( b_parity /* bond_parity_mark*/ ) {
        DrawBondParity( pDC, x1, y1, x2, y2, b_parity );
    }
    return ret;
}
/****************************************************************************/
double RoundDouble(double X)
{
    return ((X)>0.0? floor((X)+0.5):-floor(0.4999999-(X)));
}
/****************************************************************************/
int    nRound( double X )
{
    return (int)RoundDouble( X );
}
/****************************************************************************/
void roundoff_coord( double dx1, double dx2, int *new_ix1, int *new_ix2 )
{
    int new_x1;
    int new_x2;
    int diff = (int) RoundDouble(dx2-dx1);
    
    if ( dx1 >= dx2 ) {
        new_x1 = (int) RoundDouble(dx1);
        new_x2 = new_x1 + diff;
    } else {
        new_x2 = (int) RoundDouble(dx2);
        new_x1 = new_x2 - diff;
    }
    if ( new_ix1 )
        *new_ix1 = new_x1;
    if ( new_ix2 )
        *new_ix2 = new_x2;
}
/****************************************************************************/
void DrawPenColorFilledPolygon( HDC pDC, const POINT* pnt, int num )
{
    LOGPEN LogPen;
    HPEN   penNew=(HPEN)GetStockObject(BLACK_PEN);
    HPEN   penOld=(HPEN)SelectObject(pDC, penNew); /* get current pen */
    HBRUSH brushNew, brushOld;
    int    ret;
    ret = GetObject( penOld,    sizeof(LogPen), &LogPen );
    SelectObject(pDC, penOld);
    /* do not need to delete stock object penNew */
    if ( ret ) {
        if ( brushNew = CreateSolidBrush(LogPen.lopnColor) ) {
            brushOld = (HBRUSH)SelectObject( pDC, brushNew );
            Polygon( pDC, pnt, num );
            SelectObject( pDC, brushOld);
            DeleteObject( brushNew );
        }
    }
}
/****************************************************************************/
int DrawTextColorDot( HDC pDC )
{
    COLORREF clrColor = GetTextColor( pDC );
    int      width = GetSubstringWidth( pDC, 1, ".");
    HPEN     penNew   = CreatePen(PS_SOLID, 0, clrColor);
    HBRUSH   brushNew = CreateSolidBrush(clrColor);
    /*int      hfont    = GetFontHeight( pDC ); */
    POINT    pt;
    int      bSuccess = 0;
    int      nVertShift;
    int      nTextAlign = (TA_BOTTOM | TA_TOP | TA_BASELINE) & GetTextAlign( pDC );
    if ( TA_BASELINE == nTextAlign ) {
        nVertShift = 0;
    } else
    if ( TA_BOTTOM == nTextAlign ) {
        nVertShift = -GetFontDescent( pDC );
    } else
    if ( TA_TOP == nTextAlign ) {
        nVertShift = GetFontHeight( pDC ); /*GetFontAscent( pDC ); */
    } else {
        nVertShift = GetFontHeight( pDC );
    }
    if ( width % 2 ) {
        width ++;
    }

    GetCurrentPositionEx( pDC, &pt );
    if ( penNew && brushNew ) {
        HPEN   penOld   = (HPEN)  SelectObject( pDC, penNew );
        HBRUSH brushOld = (HBRUSH)SelectObject( pDC, brushNew );
        bSuccess = Ellipse( pDC, pt.x+width/2, pt.y - (3*width)/2+nVertShift, pt.x + (3*width)/2, pt.y-width/2+nVertShift );
        SelectObject( pDC, penOld );
        SelectObject( pDC, brushOld );
        if ( bSuccess ) {
            MoveToEx( pDC, pt.x+(3*width)/2, pt.y, NULL );
        }
    }
    if ( penNew )
        DeleteObject( penNew );
    if ( brushNew )
        DeleteObject( brushNew );
    return bSuccess;
}
/****************************************************************************/
int DrawBondStereo( HDC pDC, int x1, int y1, int x2, int y2, int b_stereo, int b_highlight, int bInvertBonds, COLORREF clrPen, int nPenWidth )
{
    double          lambda, mult, ax1, ay1, bond_sep, bond_width, bond_step, dx, dy, dx1, dy1;
    int             ix, ix1, iy, iy1, ix21, ix22, iy21, iy22;
    int             i, n;
    int             hfont = GetFontHeight( pDC );
    HPEN            hHighlightPen=0, hOldPen = 0;
    COLORREF        clrHighlight = CLR_MAGENTA;


    POINT           pnt[4];

    if ( b_highlight ) {
        hHighlightPen = CreatePen( PS_SOLID, nPenWidth, clrHighlight );
        if ( !hHighlightPen )
            hHighlightPen = (HPEN)GetStockObject( BLACK_PEN ); /*should not fail */
        hOldPen   = (HPEN)SelectObject( pDC, hHighlightPen );
    }
    /* the wedge (pointed) end is at (x1, y1) unless b_stereo is negative */
    if ( b_stereo < 0 ) {
        i = x1;
        x1 = x2;
        x2 = i;
        i = y1;
        y1 = y2;
        y2 = i;
        b_stereo = -b_stereo; /* 1=Up (solid triangle), 6=Down (Dashed triangle), 4=Either (zigzag triangle) */
    }
    if ( bInvertBonds ) {
        switch ( b_stereo ) {
        case STEREO_SNGL_UP:
            b_stereo = STEREO_SNGL_DOWN;
            break;
        case STEREO_SNGL_DOWN:
            b_stereo = STEREO_SNGL_UP;
            break;
        }
    }

    dx = x2 - x1;
    dy = y2 - y1;
    lambda = dx*dx + dy*dy;
    if ( lambda == 0.0 )
        return 0;
    lambda = sqrt( lambda ); /* bond length */
    
    bond_width = hfont / 6.0;
    if ( bond_width < 2.0 )
        bond_width = 2.0;  /* half-width */
    /* from x1 to x2     */
    bond_width = inchi_min( bond_width, lambda / 8.0 );
    bond_width = floor( 2.0*bond_width + 0.5 )/2;  /* half of actual bond width */
    
    bond_step = hfont / 4.0;
    bond_step = inchi_min( bond_step, lambda/6.0 );
    bond_step = inchi_max( bond_step, 3 );
    
    if ( b_stereo == 3 ) {
        bond_width *= 2.0; /* otherwise looks strange */
    }
    
    ax1 = -bond_width * dy / lambda;
    ay1 =  bond_width * dx / lambda;
    
    switch( b_stereo ) {
    
    case STEREO_SNGL_UP: /* Up */
        roundoff_coord( x2 - ax1, x2 + ax1, &ix, &ix1 );
        roundoff_coord( y2 - ay1, y2 + ay1, &iy, &iy1 );
        pnt[0].x = x1;
        pnt[0].y = y1;
        
        pnt[1].x = ix1;
        pnt[1].y = iy1;
        
        pnt[2].x = ix;
        pnt[2].y = iy;
        
        DrawPenColorFilledPolygon( pDC, pnt, 3 );
        break;
    case STEREO_DBLE_EITHER: /* cis or trans double bond */
        roundoff_coord( x2 - ax1, x2 + ax1, &ix21, &ix22 );
        roundoff_coord( y2 - ay1, y2 + ay1, &iy21, &iy22 );

        roundoff_coord( x1 - ax1, x1 + ax1, &ix, &ix1 );
        roundoff_coord( y1 - ay1, y1 + ay1, &iy, &iy1 );

        DrawLine( pDC, ix21, iy21, ix1, iy1 );
        DrawLine( pDC, ix22, iy22, ix, iy );

        break;
    
    case STEREO_SNGL_EITHER: /* either */

        n = (int)floor( lambda / bond_step );
        n = inchi_max( n, 4 );
        dx = dy = 0.0;
        ix = iy = 0;
        for ( i = 1; i <= 3*n/2; i ++ ) {
            
            mult = (double)i/(double)n/1.5;
            dx1 = (x2-x1)*mult;
            dy1 = (y2-y1)*mult;
            if ( i % 2 ) {
                dx1 += ax1*mult;
                dy1 += ay1*mult;
            } else {
                dx1 -= ax1*mult+ax1/bond_width;
                dy1 -= ay1*mult+ay1/bond_width;
            }

            roundoff_coord( dx, dx1, NULL, &ix1 );
            roundoff_coord( dy, dy1, NULL, &iy1 );
 
            /*line( win, x1 + ix, y1 + iy, x1 + ix1, y1 + iy1 ); */
            DrawLine( pDC, x1 + ix, y1 + iy, x1 + ix1, y1 + iy1 );
            dx = dx1;
            dy = dy1;
            ix = ix1;
            iy = iy1;
            /*line( win, x1 + ( int ) dx, y1 + ( int ) dy, x1 + ( int ) dx1, y1 + ( int ) dy1 ); */
        }
        break;
        
    case STEREO_SNGL_DOWN: /* Down */
        
        n = (int)floor( lambda / bond_step );
        n = inchi_max( n, 4 );
        
        for ( i = 0; i <= n; i ++ ) {
            mult = (double)i/(double)n;
            
            dx = dx1 = (x2-x1)*mult; /* x1,y1 offset */
            dy = dy1 = (y2-y1)*mult;
            
            mult = (double)(i)/(double)(n);
            
            bond_sep = ax1*mult;
            dx  += bond_sep;
            dx1 -= bond_sep+ax1/bond_width;
            
            bond_sep = ay1*mult;
            dy  += bond_sep;
            dy1 -= bond_sep+ay1/bond_width;
            
            roundoff_coord( dx, dx1, &ix, &ix1 );
            roundoff_coord( dy, dy1, &iy, &iy1 );
            
            DrawLine( pDC, x1 + ix, y1 + iy, x1 + ix1, y1 + iy1 );
        }
        break;
    }
    /* cleanup */
    if ( hOldPen )
        SelectObject( pDC, hOldPen );
    if ( hHighlightPen )
        DeleteObject( hHighlightPen );

    return 0;
}            
/****************************************************************************/
int DrawBondParity( HDC pDC, int x1, int y1, int x2, int y2, int parity_mark0 )
{
/*    int  hfont = GetFontHeight( pDC ); */
    int  xs, ys, width, height, parity_mark;
    char *p;
    COLORREF clrTextOld;

    if ( parity_mark0 < 0 ) {
        parity_mark = -parity_mark0;
        clrTextOld = SetTextColor( pDC, CLR_RED );
    } else {
        parity_mark = parity_mark0;
    }

    if ( parity_mark == BOND_MARK_ODD  )  p = "(-)"; else
    if ( parity_mark == BOND_MARK_EVEN )  p = "(+)"; else
    if ( parity_mark == BOND_MARK_UNDF )  p = "(?)"; else
    if ( parity_mark == BOND_MARK_UNKN )  p = "(u)"; else
    if ( parity_mark == BOND_MARK_ERR  )  p = "(!)"; else return 0;

    if ( abs( x2-x1 ) > 10 * abs(y2 - y1) ) {
        /* almost horizontal bond; draw parity closer (1:2) to the right end */
        if ( x2 < x1 ) {
            int tmp;
            tmp = x2;
            x2 = x1;
            x1 = tmp;
            tmp = y2;
            y2 = y1;
            y1 = tmp;
        }
        xs = 2*x1 + (4*(x2-x1))/3; /* 2/3 shift */
        ys = 2*y1 + (4*(y2-y1))/3; /* 2/3 shift */
    } else {
        xs = x1 + x2;  /* middle */
        ys = y1 + y2;  /* middle */
    }
    GetTextSize( pDC, strlen(p), p, &width, &height );
    /*
    width  = GetFontAveWidth( pDC);
    height = GetFontAscent( pDC );
    */
    xs = (xs - width  )/2;
    ys = (ys - height )/2;

    TextOut( pDC, xs, ys, p, 3 );

    if ( parity_mark0 < 0 ) {
        SetTextColor( pDC, clrTextOld );
    }

    return 0; 
}
/****************************************************************************/
int DrawBondNoStereo( HDC pDC, int x1, int y1, int x2, int y2, int b_type, int b_highlight, COLORREF clrPen, int nPenWidth )
{
    double            lambda, ax1, ax2, ay1, ay2, bond_sep;
    int             hfont = GetFontHeight( pDC );
    int             ret = 0;
    HPEN            hSolidPen=0, hDashedPen=0, hLongDashedPen=0, hHighlightPen = 0, hOldPen=0;
    char            c=0, r=0, l=0, bDashed=1, bLongDashed=0, i;
    COLORREF        clr = clrPen;

    switch( b_type ) {   /* c=center, l=left, r=right */
    case 0:
        c = 'D';
        b_highlight = 1;
        break;
    case BOND_SINGLE:    /* S=solid, D=dashed, L=long dashed(Green) */
        c = 'S';
        bDashed = 0;
        break;
    case BOND_DOUBLE:
        l = r = 'S';
        bDashed = 0;
        break;
    case BOND_TRIPLE:
        l = r = c = 'S';
        bDashed = 0;
        break;
    case BOND_TAUTOM:
        l = r = 'D';
        break;
    case BOND_ALTERN: /* 1 or 2 */
        l = 'S';
        r = 'D';
        break;
    case BOND_ALT12NS: /* 1 or 2, non-stereo */
        l = 'L';
        r = 'D';
        bLongDashed = 1;
        break;
    case BOND_ALT_13: /* 1 or 3 */
        l = 'S';
        c = r = 'D';
        break;
    case BOND_ALT_23: /* 2 or 3 */
        l = c = 'S';
        r = 'D';
        break;
    case BOND_ALT_123: /* 1 or 2 or 3 */
        c = 'S';
        l = r = 'D';
        break;
    default:
        c = 'D';
        break;
    }
    if ( bLongDashed ) {
        clr = CLR_GREEN;
    }
    /* -- debug --
    if ( c && l && r ) {
        int stop = 1;
    }
    */
    if ( b_highlight ) {
        clr = CLR_MAGENTA;
        if ( c=='S' || l == 'S' || r == 'S' ) {
            hHighlightPen = CreatePen( PS_SOLID, nPenWidth, clr );
            if ( !hHighlightPen )
                hHighlightPen = (HPEN)GetStockObject( BLACK_PEN ); /*should not fail */
            hOldPen   = (HPEN)SelectObject( pDC, hHighlightPen );
        }
    }

    if ( bDashed ) {
        hDashedPen = CreatePen( PS_DOT, 1, clr );
        if ( !hDashedPen )
            hDashedPen = (HPEN)GetStockObject( BLACK_PEN ); /*should not fail */
    }
    if ( bLongDashed ) {
        hLongDashedPen = CreatePen( PS_DASH, 1, clr );
        if ( !hLongDashedPen )
            hLongDashedPen = (HPEN)GetStockObject( BLACK_PEN ); /*should not fail */
    }

    /* draw lines between the atoms */
    if ( c == 'S' ) {
        DrawLine( pDC, x1, y1, x2, y2 );
    } else
    if ( c == 'D' ) {
        hSolidPen = (HPEN)SelectObject( pDC, hDashedPen );
        DrawLine( pDC, x1, y1, x2, y2 );
        hDashedPen = (HPEN)SelectObject( pDC, hSolidPen );
    } else
    if ( c == 'L' ) {
        hSolidPen = (HPEN)SelectObject( pDC, hLongDashedPen );
        DrawLine( pDC, x1, y1, x2, y2 );
        hLongDashedPen = (HPEN)SelectObject( pDC, hSolidPen );
    }

    /* draw lines parallel to line between bonds */
    if ( l || r ) {

        if ( b_type == BOND_DOUBLE ) {
            bond_sep = hfont / 12.0;
            if( ( x1 == x2 ) || ( y1 == y2 ) ) {
                if ( bond_sep < 1.0 ) bond_sep = 1.0;
            } else {
                if ( bond_sep < 2.0 ) bond_sep = 2.0;
            }
        } else {
            bond_sep = hfont / 6.0;
            if( ( x1 == x2 ) || ( y1 == y2 ) ) {
                if ( bond_sep < 2.0 ) bond_sep = 2.0;
            } else {
                if ( bond_sep < 4.0 ) bond_sep = 4.0;
            }
        }

        lambda = ( ( double ) ( x2 - x1 ) ) * ( x2 - x1 ) +
                 ( ( double ) ( y2 - y1 ) ) * ( y2 - y1 );
        
        if ( lambda > 0.0 )
            lambda = bond_sep / sqrt( lambda );
        else {
            ret = 1;
            goto exit_function;
        }
        for ( i = 0; i < 2; i ++, lambda = -lambda ) {
            c = i? r : l;
            ax1 = lambda * ( y1 - y2 ) + x1;
            ay1 = lambda * ( x2 - x1 ) + y1;
            ax2 = lambda * ( y1 - y2 ) + x2;
            ay2 = lambda * ( x2 - x1 ) + y2;

            if ( c == 'S' ) {
                DrawLine( pDC, ( int ) ax1, ( int ) ay1, ( int ) ax2, ( int ) ay2 );
            } else
            if ( c == 'D' ) {
                hSolidPen = (HPEN)SelectObject( pDC, hDashedPen );
                DrawLine( pDC, ( int ) ax1, ( int ) ay1, ( int ) ax2, ( int ) ay2 );
                hDashedPen = (HPEN)SelectObject( pDC, hSolidPen );
            } else
            if ( c == 'L' ) {
                hSolidPen = (HPEN)SelectObject( pDC, hLongDashedPen );
                DrawLine( pDC, x1, y1, x2, y2 );
                hLongDashedPen = (HPEN)SelectObject( pDC, hSolidPen );
            }
        }

    }

exit_function:
    /* make sure DC has its initial pen */
    if ( hSolidPen )
        SelectObject( pDC, hSolidPen );
    if ( hOldPen )
        SelectObject( pDC, hOldPen );
    /* delete all newly created pens */
    if ( hDashedPen )
        DeleteObject( hDashedPen );
    if ( hLongDashedPen )
        DeleteObject( hLongDashedPen );
    if ( hHighlightPen )
        DeleteObject( hHighlightPen );

    return ( ret );
}
/****************************************************************************/
int MoveHydrogenAtomToTheLeft( char *s, int start, int H )
{
    int  len, c, num_alpha;
    char szBuffer[16];
    char *pH = strchr( s+start, H );
    for ( pH = s+start, num_alpha = 0; (c= UCINT*pH) && c != H && c != '/'; pH ++ ) {
        num_alpha = isalpha( c );
    }
    if ( c == H && num_alpha ) { /* do not search beyond the first slash */
        for ( len = 1; pH[len] && isdigit( UCINT pH[len] ); len ++ )
            ;
        if ( len >= (int)sizeof(szBuffer) )
            return start; /* too long string */
        memcpy( szBuffer, pH, len );
        memmove( s+len, s, pH - s );
        memmove( s, szBuffer, len );
        return start+len; /* (pH-s)+i; */
    }
    return start;
}
/****************************************************************************/
int MyTextOutABC( const char *p, int iFst, int iLst, HDC pDC )
{
    ABC abc;
    POINT pt;
    if ( iFst > iLst || iFst < 0 || iLst < 0 )
        return 0;
    GetCurrentPositionEx( pDC, &pt );
    if ( GetCharABCWidths(
          pDC,   /* handle to DC */
          (int)p[iFst],  /* first character in range */
          (int)p[iFst],  /* last character in range */
          &abc   /* array of character widths */
          ) && abc.abcA < 0 ) {
        pt.x -= abc.abcA;
        MoveToEx( pDC, pt.x, pt.y, NULL );
    }
    TextOut( pDC, pt.x, pt.y, p+iFst, iLst-iFst+1 );
    if ( GetCharABCWidths(
          pDC,    /* handle to DC */
          (int)p[iLst],   /* first character in range */
          (int)p[iLst],   /* last character in range */
          &abc    /* array of character widths */
          ) && abc.abcC < 0 ) {
        GetCurrentPositionEx( pDC, &pt );
        pt.x -= abc.abcC;
        MoveToEx( pDC, pt.x, pt.y, NULL );
    }
    return 1;
}
/****************************************************************************/
int DrawColorString( HDC pDC, const char *st, int xs, int ys, int bHighlightTheAtom )
{
    int  afont = GetFontAscent( pDC );
    COLORREF clrBk = CLR_WHITE;
    int  nNumSlash = 0;
    COLORREF  OrigBkColor = CLR_WHITE;
    COLORREF  OrigTxColor = GetTextColor( pDC );
    COLORREF  NewBkColor;
    int       bWritingAtomStringColored = 0;
    /* Draw the string character by character. */
    /* For each character within the first part of the string */
    /* make a decision if it is a subscript or a superscript */
    /* The first part ends with '/'=start of the canonical number or '(' = start of the parity mark */
    /* The superscript first character is ^ or + or - or . (. initiates double shift) */
    /* The the superscript ends with not the first character and not a digit or the end of the first part */
    /* The first character of a subscript (if it is not in the superscript) is a digit */
    /* The subscript ends  a non-digit */
    const char *p=st;
    UINT uPrevTextAlign = SetTextAlign( pDC, TA_UPDATECP);
    POINT pt, pt0;
    int   i, i0, len, bSuperscript, bSubscript, bWritingAtomString, nShift, nShift2;
    int   bNewBkColor, bMoveToEx, bBypassCurrentChar;
    nShift = afont/2;
    nShift2 = nShift/2;
    if ( bHighlightTheAtom ) {
        clrBk = SetBkColor( pDC, CLR_LTPURPLE  );
    } else
    if ( strchr(st, '~') ) {
        clrBk = SetBkColor( pDC, CLR_CYAN  );
    } /* else { == for debugging ==
        clrBk = SetBkColor( pDC, CLR_LTGRAY  );
    } */
    if ( st[0] == '!' ) {
        OrigTxColor = SetTextColor( pDC, CLR_RED );
        bWritingAtomStringColored = 1;
        p ++;
    }
    bSuperscript = bSubscript = 0;
    bWritingAtomString = 1;
    MoveToEx( pDC, xs, ys, NULL );
    for ( i = i0 = 0, len = strlen( p ); i < len; i ++ ) {
        bNewBkColor = 0;
        bMoveToEx = 0;
        if ( bWritingAtomString ) { 
            GetCurrentPositionEx( pDC, &pt );
            pt0 = pt;
            if ( !p[i] || p[i] == '/' || p[i] == '(' ) {
                bWritingAtomString = 0;
                if ( bSuperscript ) {
                    while( bSuperscript > 1 ) {
                        pt.y += nShift2;
                        bSuperscript --;
                    }
                    pt.y += nShift;
                    /* MoveToEx( pDC, pt.x, pt.y, NULL ); */
                    bMoveToEx = 1;
                    bSuperscript = 0;
                }
                if ( bSubscript ) {
                    pt.y -= nShift;
                    bMoveToEx = 1;
                    /* MoveToEx( pDC, pt.x, pt.y, NULL ); */
                    bSubscript = 0;
                }
            } else {
                /* turn Superscript off */
                if ( bSuperscript == 2 && p[i] != '.' ) {
                    pt.y += nShift2;
                    bMoveToEx = 1;
                    /* MoveToEx( pDC, pt.x, pt.y, NULL ); */
                    bSuperscript -= 1;
                }
                if ( bSuperscript && !isdigit((int)(unsigned char)p[i]) && p[i] != '+' && p[i] != '-' && p[i] != '.' ) {
                    while( bSuperscript > 1 ) {
                        pt.y += nShift2;
                        bSuperscript --;
                    }
                    pt.y += nShift;
                    bMoveToEx = 1;
                    /* MoveToEx( pDC, pt.x, pt.y, NULL ); */
                    bSuperscript = 0;
                }
                /* turn Subscript off */
                if ( bSubscript && !isdigit((int)(unsigned char)p[i] ) ) {
                    pt.y -= nShift;
                    bMoveToEx = 1;
                    /* MoveToEx( pDC, pt.x, pt.y, NULL ); */
                    bSubscript = 0;
                }
                /* turn Subscript on after non-space and non-digit */
                if ( !bSuperscript && !bSubscript && ( i && p[i-1] != ' ' && !isdigit((int)(unsigned char)p[i-1]) && isdigit((int)(unsigned char)p[i])) ) {
                    pt.y += nShift;
                    bMoveToEx = 1;
                    /* MoveToEx( pDC, pt.x, pt.y, NULL ); */
                    bSubscript = 1;
                }
                /* turn Superscript on */
                if ( !bSuperscript && !bSubscript &&
                    (p[i] == '^' || p[i] == '.' ||
                    (p[i] == '+' || p[i] == '-') && (i && p[i-1]!=' ' ) ) ) {
                    pt.y -= nShift;
                    bSuperscript += 1;
                    if ( p[i] == '.' ) {
                        bSuperscript += 1; /* special case */
                        pt.y -= nShift2;
                    }
                    bMoveToEx = 1;
                    /* MoveToEx( pDC, pt.x, pt.y, NULL ); */
                    if ( p[i] == '^' ) {
                        goto output_the_string;
                        continue; /* leading superscript: isotopic mass */
                    }
                }
                if ( bSuperscript == 1 && !bSubscript && p[i] == '.' ) {
                    pt.y -= nShift2;
                    bSuperscript += 1;
                    bMoveToEx = 1;
                    /* MoveToEx( pDC, pt.x, pt.y, NULL ); */
                }

            }
        }
        if ( p[i] == '/' ) {
            nNumSlash ++;
            if ( nNumSlash == 1 ) {
                OrigBkColor = GetBkColor( pDC );
                bNewBkColor ++;
                NewBkColor = CLR_YELLOW;
                /* OrigBkColor = SetBkColor( pDC, CLR_YELLOW);*/ /* CLR_CYAN  ); */
            }
            if ( nNumSlash == 3 ) {
                bNewBkColor ++;
                NewBkColor = CLR_CYAN;
                /* SetBkColor( pDC, CLR_CYAN ); */
            }
        }
#ifdef DISPLAY_DEBUG_DATA            
        else if ( nNumSlash && p[i] == '`' ) {
            bNewBkColor ++;
            NewBkColor = CLR_LTPURPLE;
            /* SetBkColor( pDC, CLR_LTPURPLE );*/ /* DebugData */
        }
#endif
        /* output one character; special treatment for '.' */
output_the_string:
        bBypassCurrentChar = ( p[i] == '.' && bSuperscript || p[i] == '^' );
        if ( bBypassCurrentChar || bNewBkColor || bMoveToEx ) {
            /* output the accumulated string */
            int iLast = i-(bBypassCurrentChar || bMoveToEx);
            /*
            int k;
            for ( k = i0; k <= iLast; k ++ ) {
                MyTextOutABC( p, k, k, pDC );
            }
            */
            MyTextOutABC( p, i0, iLast, pDC ); /* including iLast */
            i0 = i+1-bMoveToEx;
            if ( p[i0] == '^' )
                i0 ++;
            if ( bMoveToEx ) {
                int dx = pt.x - pt0.x;
                int dy = pt.y - pt0.y;
                GetCurrentPositionEx( pDC, &pt );
                pt.x += dx;
                pt.y += dy;
                MoveToEx( pDC, pt.x, pt.y, NULL );
            }
            if ( bNewBkColor ) {
                SetBkColor( pDC, NewBkColor);
            }
            if ( p[i] == '.' && bSuperscript ) {
                if ( !DrawTextColorDot( pDC ) ) {
                    TextOut( pDC, xs, ys, p+i, 1 );
                }
            }
        }
        /*
        if ( p[i] != '.' || !bSuperscript || !DrawTextColorDot( pDC ) ) {
            TextOut( pDC, xs, ys, p+i, 1 );
        }
        */
    }
    MyTextOutABC( p, i0, i-1, pDC ); /* output the rest of the string */

    SetTextAlign( pDC, uPrevTextAlign);
    if ( nNumSlash ) {
        SetBkColor( pDC, OrigBkColor );
    }
    if ( bHighlightTheAtom || strchr(st, '~') ) {
        SetBkColor( pDC, clrBk );
    }
    if ( bWritingAtomStringColored ) {
        SetTextColor( pDC, OrigTxColor );
        bWritingAtomStringColored = 0;
    }

    return ( 1 );
}
/****************************************************************************/
int DrawPreparedString( HDC pDC, char *st1, int shift, int x, int y, int bHighlightTheAtom )
{
    DrawColorString( pDC, st1, x+shift, y-GetFontAscent( pDC )/2, bHighlightTheAtom );
    return 1;
}
/****************************************************************************/
int DrawString( HDC pDC, char *st1, int shift, int x, int y )
{
    char st[256];
    int  xs, ys, l, k;
/*    int  hfont = GetFontHeight( pDC ); */
    int  afont = GetFontAscent( pDC );
/*    COLORREF clrBk; */
/*    int  nNumSlash = 0; */
/*    COLORREF  OrigBkColor; */

    strncpy( st, st1, sizeof(st)-1 );
    st[sizeof(st)-1] = '\0';

    l = GetStringWidth( pDC, st );

    /* Default values */
    xs = x - l / 2;
    ys = y    - afont / 2;

    /* Single element */

    if( strlen( st ) == 1 ) goto draw;

    /* Changing order for right/left connection */

    if( shift  ) {
        int start = 0;
        start = MoveHydrogenAtomToTheLeft( st, start, 'T' );
        start = MoveHydrogenAtomToTheLeft( st, start, 'D' );
        start = MoveHydrogenAtomToTheLeft( st, start, 'H' );
        /* determine the position */

        if( ( strlen( st ) == 2 ) && islower( st[1] ) ) goto draw;

        k = GetStringWidth( pDC, st );
        k -= GetOneCharInStringWidth( pDC, ( st + strlen( st ) - 1 ) ) / 2;
        xs = x - k;
    } else {
        
        /* determine the position */

        k = GetOneCharInStringWidth( pDC, st );
        xs = x - k / 2;
    }

    draw:
    DrawColorString( pDC, st, xs, ys, 0 );
    return ( 1 );
}

/****************************************************************************/
void        DrawLine( HDC pDC, int x1, int y1, int x2, int y2 )
{
    MoveToEx( pDC, x1, y1, NULL );
    LineTo( pDC, x2, y2 );
}
/****************************************************************************/
int  nGetNumLegendOptions( inf_ATOM *inf_at, int num_at )
{
    int i, n, nmax=0;
    char *p;
    for ( i = 0; i < num_at; i ++ ) {
        for ( n = 0, p = inf_at[i].at_string; p = strchr( p, '/' ); p ++, n++ )
            ;
        nmax = inchi_max( nmax, n );
    }
    return nmax;
}
/****************************************************************************/
int         DrawTheInputStructure( inp_ATOM *at, INF_ATOM_DATA *inf_at_data, int num_at,
                              HDC pDC, int tx_off, int ty_off, int xoff, int yoff,
                              int width_pix, int height_pix, int bDraw, int bOrigAtom, COLORREF clrPen, int nPenWidth )
{
    int          RetVal;
    char        *NoRoom     = "Window is too small";
#ifdef INCHI_LIB
    static char  PressEnter[] = "";
#else
    static char PressEnter[] = "Press Enter to continue.";
#endif
    static char Legend[]     = "Legend:";
    char       **Str;
    int          num_str;
    inf_ATOM *inf_at = inf_at_data? inf_at_data->at : NULL;
    static char  LastString[256];
    static char  *LastStr[]     = { "Atom / Atom Id", " / Non-stereo class", " / Mobile group id", " / Mobile group class", "" };
    static char  *StrOrig[]     = {PressEnter, Legend, "Atom / Input atom number"};
    static char  *StrInfo[]     = {PressEnter, Legend, LastString};

    int          nFontHeight, nFontAveWidth, afont, i, x, y, n_opt;
    COLORREF     rgbColor;

    RetVal        = 0;
    nFontHeight   = GetFontHeight( pDC );
    nFontAveWidth = GetFontAveWidth( pDC );
    afont         = GetFontAscent( pDC );


    if ( bDraw ) {
        /* drawing */
        DrawStructure( pDC, at, inf_at_data, num_at, xoff, yoff, clrPen, nPenWidth );
        /* draw the message and legend: */
        if ( !inf_at || bOrigAtom ) {
            Str = StrOrig;
            num_str = sizeof(StrOrig)/sizeof(StrOrig[0]);
        } else {
            n_opt = nGetNumLegendOptions( inf_at, num_at );
            Str = StrInfo;
            num_str = sizeof(StrInfo)/sizeof(StrInfo[0]);
            for ( i = 0, LastString[0] = '\0'; i < n_opt && LastStr[i][0]; i ++ ) {
                strcat(LastString, LastStr[i]);
            }
        }
        rgbColor = SetBkColor( pDC, CLR_CYAN  );
        x = nFontAveWidth;
        y = nFontHeight/5;
        for ( i = 0; i < num_str; i ++ ) {
            if ( i+1 < num_str ) {
                TextOut( pDC, x+tx_off, y+ty_off, Str[i], strlen(Str[i]) );
                x += GetStringWidth( pDC, Str[i] )+2*nFontAveWidth;
                if ( i == 0 ) {
                    SetBkColor( pDC, rgbColor );
                    rgbColor = SetTextColor( pDC, CLR_BLUE /*CLR_RED*/ );
                } else
                if ( i == 1 ) {
                    SetTextColor( pDC, rgbColor );
                }
            } else {
                DrawString( pDC, Str[i], 0, x+tx_off, y+afont/2+ty_off );
            }
        }
    } else {
        rgbColor = SetTextColor( pDC, CLR_RED );
        TextOut( pDC, nFontAveWidth+tx_off, nFontHeight/5+ty_off, NoRoom, strlen(NoRoom) );
        rgbColor = SetTextColor( pDC, rgbColor );
    }
    return RetVal;
}


/****************************************************************************/
typedef struct tagTableParms {
    int    thdrHeight;
    int    thdrWidth;
    int    tcellHeight;
    int    tcellWidth;
    int    tblHeight;
    int    tblWidth;
    int    tblRows;
    int    tblCols;
    int    xtblOffs;
    int    ytblOffs;
} TBL_PARMS;
/****************************************************************************/
void CalcTblParms( HDC hMemoryDC, TBL_PARMS *tp, TBL_DRAW_PARMS *tdp,
                   int *xStructOffs, int *yStructOffs, int *xStructSize, int *yStructSize, int yoffs1)
{
    int i, j, n, w, h;
    tp->tblCols = tdp->bDrawTbl;
    tp->tblRows = 0;
    for ( i = 0; i < tp->tblCols; i ++ ) {
        (tdp->nOrientation? GetVertTextSize:GetTextSize)( hMemoryDC, strlen(tdp->ReqShownFoundTxt[i]), tdp->ReqShownFoundTxt[i], &w, &h );
        tp->thdrHeight=inchi_max(h, tp->thdrHeight); 
        tp->thdrWidth =inchi_max(w, tp->thdrWidth );

        for ( j = 0, n = 0; j < TDP_NUM_PAR; j ++ ) {
            if ( tdp->ReqShownFound[i][j] >= ' ' ) {
                GetTextSize( hMemoryDC, 1, &tdp->ReqShownFound[i][j], &w, &h );
                tp->tcellHeight=inchi_max(h, tp->tcellHeight); 
                tp->tcellWidth =inchi_max(w, tp->tcellWidth );
                n ++; /* number of types of requested or found or shown features (type: B/T, I/N, S) */
            }
        }
        tp->tblRows = inchi_max(tp->tblRows, n);
    }
    if ( tdp->nOrientation ) { /* here are tp->tblCols columns and tp->tblRows rows. */
        tp->tblHeight = tp->thdrHeight + (2*tp->tblRows+2)*tp->tcellHeight; /* empty lines above the header and around each cell */
        tp->thdrWidth = tp->tcellWidth = inchi_max( tp->tcellWidth, tp->thdrWidth );
        tp->tblWidth  = (2*tp->tblCols+1) * tp->tcellWidth; /* add empty columns around each column */
        *xStructOffs += tp->tblWidth; /* draw on the left margine */
        /* *yStructSize -= tp->tblHeight; */
        *xStructSize -= tp->tblWidth;
        tp->xtblOffs  = 0;
        tp->ytblOffs  = yoffs1;
    } else {  /* Do not believe your eyes: here are tp->tblCols rows and tp->tblRows columns. */
        tp->thdrHeight = tp->tcellHeight = inchi_max(tp->thdrHeight, tp->tcellHeight);
        tp->tblWidth   = tp->thdrWidth + (2*tp->tblRows+2)*tp->tcellWidth;
        tp->tblHeight  = (2*tp->tblCols+1)*tp->tcellHeight;
        /* draw the table on the left margine */
        *xStructOffs += tp->tblWidth;
        *xStructSize -= tp->tblWidth;
        /* *xStructSize -= tp->tblWidth; */
        tp->xtblOffs  = 0;
        tp->ytblOffs  = yoffs1;
    }
}
/****************************************************************************/
int DrawTheTable( HDC hDC, TBL_PARMS *tp, TBL_DRAW_PARMS *tdp, int x_offs, int y_offs )
{
    int i, j, ret;
    int dx = tp->tcellWidth/2;
    int dy = tp->tcellHeight/2;
    int x1, y1, x2, y2;
    /* draw frame around the table */
    ret = Rectangle( hDC, tp->xtblOffs+dx+x_offs, tp->ytblOffs+dy+y_offs, tp->xtblOffs+tp->tblWidth-dx+x_offs, tp->ytblOffs+tp->tblHeight-dy+y_offs);
    /* draw lines between labeled rows or columns */
    for ( i = 1; i < tp->tblCols; i ++ ) {
        if ( tdp->nOrientation ) {
            /* parallel to vertical columns */
            x1 = x2 = tp->xtblOffs+dx + 2 * i * tp->tcellWidth;
            y1 = tp->ytblOffs+dy;
            y2 = tp->ytblOffs+tp->tblHeight-dy;
        } else {
            /* parallel to horizontal rows */
            x1 = tp->xtblOffs+dx;
            x2 = tp->xtblOffs+tp->tblWidth-dx;
            y1 = y2 = tp->ytblOffs+dy + 2 * i * tp->tcellHeight;
        }
        DrawLine( hDC, x1+x_offs, y1+y_offs, x2+x_offs, y2+y_offs );
    }
    /* draw lines between requested/Shown/Found types */
    for ( i = 0; i < tp->tblRows; i ++ ) {
        if ( tdp->nOrientation ) {
            /* perpendicular to vertical columns */
            x1 = tp->xtblOffs+dx;
            x2 = tp->xtblOffs+tp->tblWidth-dx;
            y1 = y2 = tp->ytblOffs + tp->thdrHeight + tp->tcellHeight + 2 * i * tp->tcellHeight + dy;
        } else {
            /* perpendicular to horizontal rows */
            x1 = x2 = tp->xtblOffs + tp->thdrWidth + tp->tcellWidth + 2 * i * tp->tcellWidth + dx;
            y1 = tp->ytblOffs+dy;
            y2 = tp->ytblOffs+tp->tblHeight-dy;
        }
        DrawLine( hDC, x1+x_offs, y1+y_offs, x2+x_offs, y2+y_offs );
    }
    /* draw the text */
    for ( i = 0; i < tp->tblCols; i ++ ) {
        if ( tdp->nOrientation ) {
            /* vertical column */
            x1 = tp->xtblOffs + (2 * i + 1) * tp->tcellWidth;
            y1 = tp->ytblOffs + tp->tcellHeight;
        } else {
            /* horizontal row */
            x1 = tp->xtblOffs + tp->tcellWidth;
            y1 = tp->ytblOffs + tp->tcellHeight + 2 * i * tp->tcellHeight;
        }
        (tdp->nOrientation? TextOutVert:TextOutHoriz)( hDC, x1+x_offs, y1+y_offs, tdp->ReqShownFoundTxt[i], strlen(tdp->ReqShownFoundTxt[i]), tp->tcellWidth );
        
        for ( j = 0; j < tp->tblRows; j ++ ) {
            if ( tdp->ReqShownFound[i][j] >= ' ' ) {
                if ( tdp->nOrientation ) {
                    /* vertical column */
                    y1 = tp->ytblOffs + tp->thdrHeight + (2*j + 2) * tp->tcellHeight;
                } else {
                    /* horizontal row */
                    x1 = tp->xtblOffs + tp->thdrWidth + (2*j + 2) * tp->tcellWidth;
                }
                (tdp->nOrientation? TextOutVert:TextOutHoriz)( hDC, x1+x_offs, y1+y_offs, &tdp->ReqShownFound[i][j], 1, tp->tcellWidth );
            }
        }
    }


    return 0;
}
/****************************************************************************/
void GetStructSizes( HDC hDC, inf_ATOM *inf_at, inp_ATOM *at0, inp_ATOM *at1, int num_at, int *xoffs1, int *xoffs2, INT_DRAW_PARMS *idp)
{
    int i, j, k, num_bonds;
    double xmin, xmax, ymin, ymax;
    double x2, x, y2, y, dist;
    char  *str;
    int    len,  cur_len, half_char_width;
    int    Left_shift, Right_shift, Other_shift;
    char   cLeftChar, cRightChar;
    int    max_left_label_width_pix;
    int    max_right_label_width_pix;

    if ( idp ) {
        if ( !inf_at ) {
            idp->max_left_label_width_pix  = *xoffs1;
            idp->max_right_label_width_pix = *xoffs2;
        } else {
            idp->max_left_label_width_pix = idp->max_right_label_width_pix = 0;
        }
    } else {
        if ( !inf_at ) {
            max_left_label_width_pix  = *xoffs1;
            max_right_label_width_pix = *xoffs2;
        } else {
            max_left_label_width_pix = max_right_label_width_pix = 0;
        }
    }
    
    xmin=xmax=at0[0].x;
    ymin=ymax=at0[0].y;

    for ( num_bonds = 0, i=0; i < num_at; i ++ ) {
    
        x = at0[i].x;
        y = at0[i].y;
        Left_shift = Right_shift = Other_shift = 0;
        
        for ( j = 0; j < at0[i].valence; j ++ ) {
            k = at0[i].neighbor[j];
            x2 = at0[k].x;
            y2 = at0[k].y;
            dist = sqrt( (x-x2)*(x-x2)+(y-y2)*(y-y2) );
            if ( x < x2 - 0.2*dist )
                Left_shift ++;
            else
            if ( x > x2 + 0.1*dist )
                Right_shift ++;
            else
                Other_shift ++;
        }

        if ( inf_at ) {
            len = 0;
            str = inf_at[i].at_string;
            do {
                if ( cur_len=strcspn(str, "^") ) {
                    if ( !len ) {
                        cLeftChar = str[0];
                    }
                    cRightChar = str[cur_len-1];
                    len += GetSubstringWidth( hDC, cur_len, str);
                }
                str += cur_len+1;
            } while ( str[0] );
            inf_at[i].DrawingLabelLength = len;
            if ( Left_shift && !Right_shift ) {
                /* Atom label should be to the left from the vertex */
                half_char_width = len? GetOneCharInStringWidth( hDC, &cRightChar )/2:0;
                inf_at[i].DrawingLabelLeftShift = len - half_char_width;
                if ( idp ) {
                    idp->max_left_label_width_pix = inchi_max( idp->max_left_label_width_pix, inf_at[i].DrawingLabelLeftShift);
                    idp->max_right_label_width_pix = inchi_max(idp->max_right_label_width_pix, half_char_width);
                } else {
                    max_left_label_width_pix  = inchi_max(max_left_label_width_pix, inf_at[i].DrawingLabelLeftShift);
                    max_right_label_width_pix = inchi_max(max_right_label_width_pix, half_char_width);
                }
                /* convert NH2 to H2N, etc. */
                len = 0;
                str = inf_at[i].at_string;
                len = MoveHydrogenAtomToTheLeft( str, len, 'T' );
                len = MoveHydrogenAtomToTheLeft( str, len, 'D' );
                len = MoveHydrogenAtomToTheLeft( str, len, 'H' );
            } else {
                /* Atom label should be to the right from the vertex */
                half_char_width = len? GetOneCharInStringWidth( hDC, &cLeftChar )/2:0;
                inf_at[i].DrawingLabelLeftShift = half_char_width;
                if ( idp ) {
                    idp->max_left_label_width_pix = inchi_max( idp->max_left_label_width_pix, half_char_width);
                    idp->max_right_label_width_pix = inchi_max( idp->max_right_label_width_pix,
                                                          inf_at[i].DrawingLabelLength
                                                        - inf_at[i].DrawingLabelLeftShift);
                } else {
                    max_left_label_width_pix  = inchi_max( max_left_label_width_pix, half_char_width);
                    max_right_label_width_pix = inchi_max( max_right_label_width_pix,
                                                          inf_at[i].DrawingLabelLength
                                                        - inf_at[i].DrawingLabelLeftShift);
                }
            }
        } else {
            at1[i].bDrawingLabelLeftShift = ( Left_shift && !Right_shift );
        }

        xmin = inchi_min( xmin, x );
        xmax = inchi_max( xmax, x );
        ymin = inchi_min( ymin, y );
        ymax = inchi_max( ymax, y );
    }
    if ( idp ) {
        idp->xmin = xmin;
        idp->xmax = xmax;
        idp->ymin = ymin;
        idp->ymax = ymax;
        if ( inf_at ) {
            *xoffs1 = idp->max_left_label_width_pix;
            *xoffs2 = idp->max_right_label_width_pix;
        }
    } else
    if ( inf_at ) {
        *xoffs1 = max_left_label_width_pix;
        *xoffs2 = max_right_label_width_pix;
    }

}
/****************************************************************************/
void ResizeAtomForDrawing( inf_ATOM *inf_at, inp_ATOM *at0, inp_ATOM *at1, int num_at,
                          INT_DRAW_PARMS *idp, int width, int height, int nFontWidth, int *xoffs1, int *xoffs2,
                          int *draw_width, int *draw_height, int *xdraw_offs, int *ydraw_offs )
{ 
    int    i;
    double xmin = idp->xmin;
    double xmax = idp->xmax;
    double ymin = idp->ymin;
    double ymax = idp->ymax;
    double dx = 0.0, dy = 0.0, new_dx;
    double coeff=0.0, xShift=0.0, yShift=0.0;
    double coeffx = 0.0, coeffy = 0.0, new_coeffx;


    if ( xmax > xmin || ymax > ymin ) {

        dx = xmax-xmin;
        dy = ymax-ymin;

        if ( width > 0 && height > 0 ) {
            coeffx = dx > 0.0? (double)width/dx  : 0.0;
            coeffy = dy > 0.0? (double)height/dy : 0.0;

            if ( coeffx > 0.0 && coeffy > 0.0 )
                coeff = inchi_min( coeffx, coeffy );
            else
                coeff = inchi_max( coeffx, coeffy );
        }
        if ( coeffx == 0.0 ) {
            xShift = width/2.0;
        }
        if ( coeffy == 0.0 ) {
            yShift = height/2.0;
        }

    } else {
        coeff = 0.0;
        xShift = width/2.0;
        yShift = height/2.0;
    }
   
    
    /* set screen coordinates for drawing */
    for ( i = 0; i < num_at; i ++ ) {
        at1[i].y = (ymax - at0[i].y)*coeff+yShift; /* screen y axis is directed down */
        at1[i].x = (at0[i].x-xmin)*coeff+xShift;
    }
    /* horizontal screen coordinates rescale if x dimension determines struct. size */
    if ( coeffx > 0.0 && coeffy > 0.0 && inf_at ) { 
        double new_xmin =  1.0e32;
        double new_xmax = -1.0e32;
        double dx1, dx2;
        int    nPass=0;
        int    new_width = width + idp->max_left_label_width_pix + idp->max_right_label_width_pix;
        for ( i = 0; i < num_at; i ++ ) {
            dx1 = at1[i].x - inf_at[i].DrawingLabelLeftShift;
            dx2 = dx1 + inf_at[i].DrawingLabelLength;
            new_xmin = inchi_min( new_xmin, dx1 );
            new_xmax = inchi_max( new_xmax, dx2 );
        }
        new_dx = new_xmax - new_xmin;
        if ( coeffx > coeffy && new_dx < (double)new_width )
            goto done;
        if ( new_dx < (double)new_width && new_dx > (double)(new_width-2*nFontWidth) )
            goto done;
again:
        new_dx += nFontWidth; /* precaution */
        new_coeffx = coeffx + ((double)new_width-new_dx)/dx;
        if ( coeffy > 0.0 )
            coeff = inchi_min(new_coeffx, coeffy);
        else
            coeff = new_coeffx;
        
        new_xmin =  1.0e32;
        new_xmax = -1.0e32;
        for ( i = 0; i < num_at; i ++ ) {
            dy = at1[i].y = (ymax - at0[i].y)*coeff+yShift; /* screen y axis is directed down */
            dx = at1[i].x = (at0[i].x-xmin)*coeff+xShift;
            dx1 = dx  - inf_at[i].DrawingLabelLeftShift;
            dx2 = dx1 + inf_at[i].DrawingLabelLength;
            new_xmin = inchi_min( new_xmin, dx1 );
            new_xmax = inchi_max( new_xmax, dx2 );
        }
        new_dx = new_xmax - new_xmin;
        if ( new_dx > new_width && nPass++ < 3 ) {
            coeffx = new_coeffx;
            goto again;
        }
done:
        *xoffs1      = -nRound( new_xmin );
        *xoffs2      =  nRound( new_xmax - (xmax-xmin)*coeff );
        *draw_width  =  nRound( (xmax-xmin)*coeff ); /* nRound( new_dx ); */
        *draw_height =  nRound( (ymax-ymin)*coeff );
    } else {
        *draw_width = nRound( (xmax-xmin)*coeff );
        *draw_height = nRound( (ymax-ymin)*coeff );
    }
    *xdraw_offs = nRound(xShift);
    *ydraw_offs = nRound(yShift);
}

/******************************************************************************/
void InpStructureMarkEquComponents( MY_WINDOW_DATA *pWinData, AT_NUMB nNewEquLabel,
                                   inp_ATOM *at0, inp_ATOM *at1, inf_ATOM *inf_at, int num_at  )
{ 
    int bHighlight = 0;
    AT_NUMB *nEquLabels   = pWinData->nEquLabels;
    /* highlight equivalent components */
    int i, neigh, j, ni, nj, nh;
    if ( nNewEquLabel ) {
        for ( i = 0; i < num_at; i ++ ) {
            ni = (int)at0[i].orig_at_number-1;
            if ( 0 <= ni && ni < num_at ) {
                if ( nNewEquLabel == nEquLabels[ni] ||
                     1==at0[i].el_number && 1==at0[i].chem_bonds_valence &&
                     0<=(nh=(int)at0[at0[i].neighbor[0]].orig_at_number-1) &&
                     nh < num_at && nNewEquLabel == nEquLabels[nh]  ) {
                    inf_at[i].cHighlightTheAtom = 1; /* highlight the atom */
                    bHighlight |= 1;
                    for ( j = 0; j < at0[i].valence; j ++ ) {
                        neigh = (int)at0[i].neighbor[j];
                        if ( neigh < num_at &&
                             0 <= (nj = (int)at0[neigh].orig_at_number-1) &&
                             nj < num_at &&
                             ( 
                                /* highlighted atom */
                              ( nNewEquLabel == nEquLabels[nj] ) ||
                                 /* terminal H */
                              ( 1==at0[neigh].el_number && 1 == at0[neigh].chem_bonds_valence)
                             )
                            ) {
                            at0[i].bond_type[j] |= BOND_MARK_HIGHLIGHT; /* highlight the bond */
                            at1[i].bond_type[j] |= BOND_MARK_HIGHLIGHT; /* highlight the bond */
                        }
                    }
                } else
                if ( inf_at[i].cHighlightTheAtom ) {
                    inf_at[i].cHighlightTheAtom = 0;
                    for ( j = 0; j < at0[i].valence; j ++ ) {
                            at0[i].bond_type[j] &= ~BOND_MARK_HIGHLIGHT; /* remove highlight from the bond */
                            at1[i].bond_type[j] &= ~BOND_MARK_HIGHLIGHT; /* remove highlight from the bond */
                    }

                }
            }
        }
    } else {
        for ( i = 0; i < num_at; i ++ ) {
            ni = (int)at0[i].orig_at_number-1;
            if ( 0 <= ni && ni < num_at ) {
                if ( inf_at[i].cHighlightTheAtom ) {
                    inf_at[i].cHighlightTheAtom = 0;
                    for ( j = 0; j < at0[i].valence; j ++ ) {
                            at0[i].bond_type[j] &= ~BOND_MARK_HIGHLIGHT; /* remove highlight from the bond */
                            at1[i].bond_type[j] &= ~BOND_MARK_HIGHLIGHT; /* remove highlight from the bond */
                    }

                }

            }
        }
    }
    if ( !bHighlight ) {
        nNewEquLabel = 0;
    }
    pWinData->nCurEquLabel = nNewEquLabel;
    pWinData->bHighlight   = bHighlight;
}
/****************************************************************************/
int CreateInputStructPicture( HDC hDC, MY_WINDOW_DATA *pWinData, RECT *rc, int bPrint, AT_NUMB nNewEquLabel )
{
    int              ErrCode = 1, Res, width=0, height=0, yoffs0=0, xoffs1=0, yoffs1=0, xoffs2, yoffs2;
    int              xDim, yDim, w, h, xs, ys;

    HDC hMemoryDC     = NULL;
    HBITMAP hBitmap  = NULL, hOldBitmap=NULL;
    HPEN    Pen      = 0, OldPen = 0;
    LOGFONT MyLogFont;
    HFONT    Font     = 0, OldFont = 0;
    char   *FaceName = FONT_NAME;
    int     bDrawTbl = 0;
    int     bStereoFlags = 0;

    int win_top         =  rc->top;
    int win_left        =  rc->left;

    int win_width       =  rc->right - rc->left;
    int win_height      =  rc->bottom - rc->top;

    int bm_top         =  rc->top;
    int bm_left        =  rc->left;

    int bm_width       =  rc->right - rc->left;
    int bm_height      =  rc->bottom - rc->top;

    TBL_PARMS       tp;

    int             num_at    =  0;
    int             bOrigAtom =  0;
    int             nFontSize =  10;
    inp_ATOM       *at0       =  NULL;
    inp_ATOM       *at1       =  NULL;
    inf_ATOM       *inf_at    =  NULL;
    INT_DRAW_PARMS *idp       =  NULL;
    TBL_DRAW_PARMS *tdp       =  NULL;
    INT_DRAW_PARMS  idp_print;

    /* structure + headers rect: offsets, width, height */
    int    xStructOffs=0, yStructOffs=0, xStructSize, yStructSize;

    int nFontHeight=0;
    int nFontWidth=0;
    int nPenWidth = 1;
    COLORREF clrPen = CLR_BLUE;
    int afont;

    /*bPrint = 1;*/ /* test */

    if ( pWinData ) {
        num_at    =  pWinData->num_at;
        bOrigAtom =  pWinData->bOrigAtom;
        nFontSize =  pWinData->nFontSize;
        at0       =  pWinData->at0;
        at1       =  pWinData->at1;
        inf_at    =  pWinData->inf_at_data.at;
        idp       = &pWinData->idp;
        tdp       = &pWinData->tdp;
        bStereoFlags = pWinData->inf_at_data.StereoFlags;
        if ( bPrint ) {
            idp = &idp_print;
            memset( idp, 0, sizeof(idp[0]) );
        }

        if ( pWinData->nCurEquLabel != nNewEquLabel &&
             pWinData->nEquLabels &&
             nNewEquLabel <= pWinData->nNumEquSets && at0 && at1 && inf_at && num_at ) {
            InpStructureMarkEquComponents( pWinData, nNewEquLabel, at0, at1, inf_at, num_at );
        }
    }


    xDim = xStructSize = win_width;
    yDim = yStructSize = win_height;
    if ( !bPrint ) {
        /* create bitmap: drawing on a bitmap reduces screen flicker */
        if ( !(hMemoryDC = CreateCompatibleDC(hDC)) ||
             !(hBitmap   = CreateCompatibleBitmap(hDC, xDim, yDim )) ||
             !(hOldBitmap = (HBITMAP) SelectObject(hMemoryDC, hBitmap)) ||
             !PatBlt( hMemoryDC, 0, 0, xDim, yDim, PATCOPY )) {
            ErrCode = 0;
            goto _end;
        }
        bm_top  = 0;
        bm_left = 0;
        bm_height = win_height - win_top;
        bm_width  = win_width  - win_left;
    } else {
        hMemoryDC = hDC;
    }

    if ( pWinData ) {
        /* create drawing tools: font */
        memset( &MyLogFont, 0, sizeof( LOGFONT ) );
        if ( nFontSize < 0 ) {
            int iLogPixsY = GetDeviceCaps(hDC, LOGPIXELSY);
            nFontSize = -MulDiv(iLogPixsY, -nFontSize, 72);
        }
        nPenWidth = bPrint? inchi_max(abs(nFontSize)/10,1):1;
        clrPen    = bPrint? CLR_BLACK : CLR_BLUE;
        MyLogFont.lfHeight = nFontSize;
        MyLogFont.lfWeight = FW_NORMAL;
        /* MyLogFont.lfItalic = 1; */ /* test MyTextOutABC() */
        strncpy( MyLogFont.lfFaceName, FaceName, LF_FACESIZE );
        Font = CreateFontIndirect( &MyLogFont ); /* black */

        /* create drawing tools: pen */
        Pen = CreatePen( PS_SOLID, nPenWidth, clrPen );

        /* select drawing tools into the bitmap */
        OldFont = (HFONT)SelectObject( hMemoryDC, Font );
        OldPen  = (HPEN) SelectObject( hMemoryDC, Pen );

        /* find sizes */

        nFontHeight = GetFontHeight( hMemoryDC );
        nFontWidth  = GetFontAveWidth( hMemoryDC );
        afont       = GetFontAscent( hMemoryDC );
    
        /* offsets within the (0, 0, xStructSize, yStructSize) rectangle */
        xoffs1 = xoffs2 = 16*nFontWidth;     /* define structure atom coordinate margins here */
        yoffs0 = (bPrint && pWinData->szTitle && pWinData->szTitle[0] )? (3*nFontHeight)/2 : 0; /* 1.5 or 0 lines at the top */
        yoffs1 = (5*nFontHeight)/2; /* 2.5 lines at the top */
        yoffs2 = (5*nFontHeight)/2; /* 2.5 lines at the bottom */

        /***********************************************/
        /*         Calculate structure size            */
        /***********************************************/
        if ( idp->bInit ) {
            /* structure sizes are known */
            xoffs1 = idp->max_left_label_width_pix;
            xoffs2 = idp->max_right_label_width_pix;
        } else {
            GetStructSizes( hMemoryDC, inf_at, at0, at1, num_at, &xoffs1, &xoffs2, idp);
            idp->bInit = 1;
        }
    
        /***********************************************/
        /* Calculate requested/Shown/Found table sizes */
        /***********************************************/
    
        memset( &tp, 0, sizeof(tp) );
#ifdef INCHI_LIB
        bDrawTbl = 0;
#else
        bDrawTbl = tdp && tdp->bDrawTbl;
        /*bDrawTbl = 0;*/
#endif
        if ( bDrawTbl ) {
            double dx = idp->xmax - idp->xmin;
            double dy = idp->ymax - idp->ymin;
            int nOrientation_tmp = tdp->nOrientation;
            if ( dx > 0.0 && dy > 0.0 ) {
                int xStructOffs_tmp0=xStructOffs, yStructOffs_tmp0=yStructOffs;
                int xStructSize_tmp0=xStructSize, yStructSize_tmp0=yStructSize;
                int twidth0;
                /*
                int xStructOffs_tmp1=xStructOffs, yStructOffs_tmp1=yStructOffs;
                int xStructSize_tmp1=xStructSize, yStructSize_tmp1=yStructSize;
                int twidth1;
                */
                tdp->nOrientation = 0;
                CalcTblParms( hMemoryDC, &tp, tdp,
                              &xStructOffs_tmp0, &yStructOffs_tmp0, &xStructSize_tmp0, &yStructSize_tmp0, yoffs0+yoffs1);
                twidth0 = tp.tblWidth;
                xStructSize_tmp0 -= xoffs1 + xoffs2;
                yStructSize_tmp0 -= yoffs1 + yoffs2;
                /*
                memset( &tp, 0, sizeof(tp) );
                tdp->nOrientation = 0;
                CalcTblParms( hMemoryDC, &tp, tdp,
                              &xStructOffs_tmp1, &yStructOffs_tmp1, &xStructSize_tmp1, &yStructSize_tmp1, yoffs1);
                twidth1 = tp.tblWidth;
                xStructSize_tmp1 -= xoffs1 + xoffs2;
                yStructSize_tmp1 -= yoffs1 + yoffs2;
                */
                if ( xStructSize_tmp0 > 0 && yStructSize_tmp0 > 0 ) {
                    nOrientation_tmp = ( (double)yStructSize_tmp0/(double)xStructSize_tmp0 > dy / dx );
                }
            }
            tdp->nOrientation = nOrientation_tmp;
            memset( &tp, 0, sizeof(tp) );
            CalcTblParms( hMemoryDC, &tp, tdp,
                          &xStructOffs, &yStructOffs, &xStructSize, &yStructSize, yoffs0+yoffs1);
        } else {
            xStructSize -= 2*nFontWidth;  /* drawing area sizes */
            yStructSize -= nFontHeight;
            xStructOffs += nFontWidth;    /* drawing area offsets */
            yStructOffs += nFontHeight/2;
        }

        width  = xStructSize - xoffs1 - xoffs2; /* drawing structre area sizes */
        height = yStructSize - yoffs0 - yoffs1 - yoffs2;
        ResizeAtomForDrawing( inf_at, at0, at1, num_at, idp, width, height, nFontWidth, &xoffs1, &xoffs2, &w, &h, &xs, &ys );
    
        /* At this point xStructOffs = the left margin for the drawing */
    
        if ( 2*xs+xoffs1+xoffs2+w < win_width - xStructOffs ) { /*compare 2*x-coordinate of the center */
            xStructOffs += (win_width - (xStructOffs+xoffs1+xoffs2+w))/2-xs;
        }
    
    }
    
    /* draw */
    if ( bPrint ) {
        ; /*PatBlt( hMemoryDC, 0, 0, xDim+win_left, yDim+win_top, PATCOPY ); */
    } else {
        PatBlt( hMemoryDC, 0, 0, xDim, yDim, PATCOPY );
    }

    if ( pWinData ) {
          char str[128]="";

        /* exact rectangle around the structure drawing */
        /*Rectangle( hMemoryDC, xStructOffs+xs-1, yStructOffs+ys+yoffs1-afont-1, xStructOffs+xs+xoffs1+xoffs2+w+1, yStructOffs+ys+yoffs1+nFontHeight+h+1); */
        Res = DrawTheInputStructure( at1, &pWinData->inf_at_data, num_at, hMemoryDC,
                                     bm_left,                        /* text output offsets */
                                     bm_top+yoffs0,
                                     bm_left + xoffs1+xStructOffs,   /* structure offsets */
                                     bm_top  + +yoffs0+yoffs1+yStructOffs,
                                     xDim-xoffs1,                    /* structure width */
                                     yDim-yoffs0-yoffs1,             /* structure height */
                                     ( width >= 0 && height >= 0 ), bOrigAtom, clrPen, nPenWidth );
        if( Res == -1 ){
            ErrCode = 0;
            goto _end;
        }
        if ( bDrawTbl || bStereoFlags || pWinData->inf_at_data.szRemovedProtons[0] ) {
            /*
            if ( bDrawTbl ) {
                char *str = "Abbreviations: Tautomeric, Isotopic, Stereo";
                DrawTheTable( hMemoryDC, &tp, tdp, bPrint?win_left:0, bPrint?win_top:0 );
                TextOut( hMemoryDC, nFontWidth+(bPrint?win_left:0), win_height - nFontHeight+(bPrint?win_top:0), str, strlen(str) ); 
            }
            */
            if ( bStereoFlags ) {
                switch ( bStereoFlags & INF_STEREO_ABS_REL_RAC ) {
                case INF_STEREO_ABS:
                    strcat( str, "Absolute" );
                    break;
                case INF_STEREO_REL:
                    strcat( str, "Relative" );
                    break;
                case INF_STEREO_RAC:
                    strcat( str, "Racemic mixture" );
                    break;
                }
                if ( str[0] ) {
                    strcat( str, " stereo" );
                    switch( bStereoFlags & INF_STEREO_NORM_INV ) {
                    case INF_STEREO_NORM:
                        strcat( str, " (normal)" );
                        break;
                    case INF_STEREO_INV:
                        strcat( str, " (inverted)" );
                        break;
                    case INF_STEREO_NORM_INV:
                        strcat( str, " (normal and inverted)" );
                        break;
                    }
                }
            }
            if ( bDrawTbl ) {
                int bTaut=0, bIso=0, bSter=0;
                if ( str[0] ) {
                    strcat( str, ";  ");
                }
                bTaut = (NULL != memchr(tdp->ReqShownFound[ilSHOWN], 'T', TDP_NUM_PAR));
                bIso  = (NULL != memchr(tdp->ReqShownFound[ilSHOWN], 'I', TDP_NUM_PAR));
                bSter = (NULL != memchr(tdp->ReqShownFound[ilSHOWN], 'S', TDP_NUM_PAR)) ||
                        (NULL != memchr(tdp->ReqShownFound[ilSHOWN], 's', TDP_NUM_PAR));
                strcat( str, "Abbreviation" );
                if ( bTaut+bIso+bSter > 1 ) {
                    strcat( str, "s:");
                } else {
                    strcat( str, ":" );
                }
                if ( bTaut ) strcat( str, " Mobile H" );
                if ( bIso  ) strcat( str, " Isotopic" );
                if ( bSter ) strcat( str, " Stereo" );
                DrawTheTable( hMemoryDC, &tp, tdp, bm_left, bm_top);
            }
            if ( pWinData->inf_at_data.szRemovedProtons[0] ) {
                if ( str[0] ) strcat( str, ";    ");
                strcat( str, pWinData->inf_at_data.szRemovedProtons );
            }
            /*TextOut( hMemoryDC, nFontWidth+bm_left, bm_height - nFontHeight +bm_top, str, strlen(str) );*/
            /*DrawColorString( hMemoryDC, str, nFontWidth+bm_left, bm_height - nFontHeight +bm_top, 0 );*/
        }
        if ( pWinData->bHighlight ) {
            /* draw highlighted (identical) components description */
            char     *p1             = " Highlighted ";
            char     *p2             = " components are identical";
            int       x              = bm_left + nFontWidth;
            int       y              = bm_top  + bm_height - nFontHeight;
            COLORREF  clrBk;
            UINT      uPrevTextAlign;
            POINT     pt;
            /* save current position */
            GetCurrentPositionEx( hMemoryDC, &pt );
            /* move to the starting point */
            MoveToEx( hMemoryDC, x, y, NULL );
            /* set text aligh that do not require coordinates in TextOut() */
            uPrevTextAlign = SetTextAlign( hMemoryDC, TA_UPDATECP);
            /* set highlighted background color */
            clrBk = SetBkColor( hMemoryDC, CLR_LTPURPLE  );
            /* output the 1st word */
            TextOut( hMemoryDC, 0, 0, p1, strlen(p1) );
            /* restore text background */
            SetBkColor( hMemoryDC, clrBk );
            /* output the rest of the text as normal text */
            TextOut( hMemoryDC, 0, 0, p2, strlen(p2) );

            if ( str[0] ) {
                POINT     pt2;
                TextOut( hMemoryDC, 0, 0, ";", 1 );
                GetCurrentPositionEx( hMemoryDC, &pt2 );
                DrawColorString( hMemoryDC, str, pt2.x+2*nFontWidth, pt2.y, 0 );
            }
            /* restore text align */
            SetTextAlign( hMemoryDC, uPrevTextAlign);
            /* restore current position */
            MoveToEx( hMemoryDC, pt.x, pt.y, NULL );
        } else
        if ( str[0] ) {
            DrawColorString( hMemoryDC, str, nFontWidth+bm_left, bm_height - nFontHeight +bm_top, 0 );
        }
        if ( yoffs0 > 0 && bPrint && pWinData->szTitle && pWinData->szTitle[0] ) {
            /* print window title */
            char     *p1             = pWinData->szTitle;
            int       x              = bm_left + 3*nFontWidth;
            int       y              = bm_top  + yoffs0 - (3*nFontHeight)/2;
            UINT      uPrevTextAlign;
            POINT     pt;
            /* save current position */
            GetCurrentPositionEx( hMemoryDC, &pt );
            /* move to the starting point */
            MoveToEx( hMemoryDC, x, y, NULL );
            /* set text aligh that do not require coordinates in TextOut() */
            uPrevTextAlign = SetTextAlign( hMemoryDC, TA_UPDATECP);
            /* output the text */
            TextOut( hMemoryDC, 0, 0, p1, strlen(p1) );
            /* restore text align */
            SetTextAlign( hMemoryDC, uPrevTextAlign);
            /* restore current position */
            MoveToEx( hMemoryDC, pt.x, pt.y, NULL );
        }
    }

    if ( !bPrint ) {
        /* copy bitmap onto the window */
        ErrCode = BitBlt(
            hDC,        /* handle to the destination device context  */
            win_left,    /* x-coordinate of destination rectangle's upper-left corner */
            win_top,    /* y-coordinate of destination rectangle's upper-left corner */
            xDim,        /* width of destination rectangle  */
            yDim,       /* height of destination rectangle  */
            hMemoryDC,    /* handle to source device context  */
            bm_left,    /* x-coordinate of source rectangle's upper-left corner   */
            bm_top,        /* y-coordinate of source rectangle's upper-left corner */
            SRCCOPY     /* raster operation code  */
        );
    }

    if ( pWinData ) {

        /* remove drawing tools */
        if( Pen ){
            SelectObject( hMemoryDC, OldPen );
            DeleteObject( Pen );
        }
        if( Font ){
            SelectObject( hMemoryDC, OldFont );
            DeleteObject( Font );
        }

    }
_end:
    if ( !bPrint ) {
        if( hBitmap ) {
            if ( hOldBitmap )
                SelectObject(hMemoryDC, hOldBitmap);
            DeleteObject( hBitmap );
        }
        if( hMemoryDC && hMemoryDC != hDC )
            DeleteDC(hMemoryDC);
    }
    return ErrCode;
}

/*********************************************************************

  FUNCTION: WndProcDisplayCanonStructure(HWND, unsigned, WORD, LONG)

  PURPOSE:  Processes messages for the main window.

  MESSAGES:

    WM_COMMAND - process the application menu
    WM_PAINT - Paint the main window
    WM_DESTROY - post a quit message and return
    WM_DISPLAYCHANGE - message sent to Plug & Play systems when the display changes
    WM_RBUTTONDOWN - Right mouse click -- put up context menu here if appropriate
    WM_NCRBUTTONUP - User has clicked the right button on the application's system menu

********************************************************************/
LRESULT CALLBACK WndProcDisplayInputStructure(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    int wmId, wmEvent;
    PAINTSTRUCT ps;
    HDC hdc;
    MY_WINDOW_DATA *pWinData;
    RECT rc;

#define IS_WIN95 0

    switch (message) { 

        case WM_COMMAND:
            wmId    = LOWORD(wParam); /* Remember, these are... */
            wmEvent = HIWORD(wParam); /* ...different for Win32! */
            break;

        case WM_CLOSE:
            pWinData = (MY_WINDOW_DATA *)GetWindowLong( hWnd, GWL_USERDATA );
            pWinData->bEsc = 1;
            goto close_window;

        case WM_RBUTTONUP: /* RightClick in the window client area */
        case WM_LBUTTONUP: /* LeftClick in the window client area */
                /* stop the timer */
                pWinData = (MY_WINDOW_DATA *)GetWindowLong( hWnd, GWL_USERDATA );
                if ( pWinData->nTimerId ) {
                    KillTimer( hWnd, pWinData->nTimerId );
                    pWinData->nTimerId = 0;
                    pWinData->bUserIntervened = 1;
                }
                /*InvalidateRect( hWnd, NULL, 0 ); */
                break;
            
        case WM_CHAR:
            
                pWinData = (MY_WINDOW_DATA *)GetWindowLong( hWnd, GWL_USERDATA );
                if ( pWinData->nTimerId ) {
                    KillTimer( hWnd, pWinData->nTimerId );
                    pWinData->nTimerId = 0;
                    pWinData->bUserIntervened = 1;
                }
                if ( wParam == '\r' || wParam == 27 ) {
                    pWinData->bEsc = (wParam == 27);
                    goto close_window;
                }
            
                break;

        case WM_ERASEBKGND:
                return TRUE; /* to prevent flicker do not let Windows erase background */

        case WM_SIZE:
        case WM_MOVE:
                pWinData = (MY_WINDOW_DATA *)GetWindowLong( hWnd, GWL_USERDATA );
                if ( pWinData->nTimerId ) {
                    KillTimer( hWnd, pWinData->nTimerId );
                    pWinData->nTimerId = 0;
                    pWinData->bUserIntervened = 1;
                }
                break;

        case WM_DISPLAYCHANGE: /* Only comes through on plug'n'play systems */
            {
                SIZE szScreen;
                BOOL fChanged = (BOOL)wParam;

                szScreen.cx = LOWORD(lParam);
                szScreen.cy = HIWORD(lParam);
                
                if (fChanged) {
                    /* The display 'has' changed. szScreen reflects the */
                    /* new size. */
                    ; /*MessageBox (GetFocus(), "Display Changed", szWindowClassName, 0); */
                } else {
                    /* The display 'is' changing. szScreen reflects the */
                    /* original size. */
                    MessageBeep(0);
                }
            }
        break;

        case WM_TIMER:
            pWinData = (MY_WINDOW_DATA *)GetWindowLong( hWnd, GWL_USERDATA );
            if ( wParam == pWinData->nTimerId ) {
                KillTimer( hWnd, pWinData->nTimerId );
                pWinData->nTimerId = 0;
                goto close_window;
            }
            break;
/*
        case WM_SHOWWINDOW:
            break;
*/
        case WM_PAINT:
            pWinData = (MY_WINDOW_DATA *)GetWindowLong( hWnd, GWL_USERDATA );
            GetClientRect( hWnd, &rc );
            hdc = BeginPaint (hWnd, &ps);
            /* the drawing code is here */
            /*
            rc.top = 30;
            rc.left= 50;
            */
            CreateInputStructPicture( hdc, pWinData, &rc, 0, pWinData->nNewEquLabel );
            EndPaint (hWnd, &ps);
            /* start the timer if requested */
            if ( !pWinData->nTimerId && !pWinData->bUserIntervened && pWinData->ulDisplTime ) {
                pWinData->nTimerId = SetTimer(
                                          hWnd,                   /* handle to window */
                                          MY_TIMER_ID,            /* timer identifier */
                                          pWinData->ulDisplTime,  /* time-out value */
                                          NULL                    /* ptr to the timer procedure */
                                        );
            }
            break;

        case WM_DESTROY:
            /* Tell WinHelp we don't need it any more... */
            /*    WinHelp (hWnd, APPNAME".HLP", HELP_QUIT,(DWORD)0); */
            /*PostQuitMessage(0); */
            return (DefWindowProc(hWnd, message, wParam, lParam));
            break;

        default:
            return (DefWindowProc(hWnd, message, wParam, lParam));
    }
    goto exit_function;

close_window:

    pWinData = (MY_WINDOW_DATA *)GetWindowLong( hWnd, GWL_USERDATA );
    GetWindowRect( hWnd, &pWinData->rc );
    ReallySetForegroundWindow(GetConsoleHwnd());
    DestroyWindow( hWnd );


exit_function:

    return (0);
}
/**********************************************************************************************/
void FreeWinData( MY_WINDOW_DATA* pWinData )
{
    if ( pWinData ) {
        if ( pWinData->at0 ) {
            inchi_free( pWinData->at0 );
            pWinData->at0 = NULL;
        }
        if ( pWinData->at1 ) {
            inchi_free( pWinData->at1 );
            pWinData->at1 = NULL;
        }
        FreeInfoAtomData( &pWinData->inf_at_data );

        if ( pWinData->nEquLabels ) {
            inchi_free( pWinData->nEquLabels );
            pWinData->nEquLabels = NULL;
            pWinData->nNumEquSets = 0;
        }

        if ( pWinData->szTitle ) {
            inchi_free(pWinData->szTitle);
            pWinData->szTitle = NULL;
        }
    }
}
#ifndef INCHI_LIB
/**********************************************************************************************/
int DisplayInputStructure( char *szOutputString, inp_ATOM  *at, INF_ATOM_DATA *inf_at_data, int num_at, DRAW_PARMS *dp )
{
#define IS_WIN95 0
    HWND hWnd = NULL;
    WNDCLASS  wc;
    HINSTANCE hInstance = 0;
    RECT      rc, rc2, rc3;
    MSG       msg;
    MY_WINDOW_DATA WinData;
    int       ret, ret2, ret3, bRectVisible;
    HDC       hDesktopDC;
    inf_ATOM  *inf_at = inf_at_data? inf_at_data->at : NULL;
    /*WINDOWPLACEMENT wndpl = {sizeof(WINDOWPLACEMENT),}; */ /*set wndpl.length */
    /* get console application window handle */
    HWND hConsoleWnd = GetConsoleHwnd();
    HWND hDesktopWnd = GetDesktopWindow();
    int  bSetForeground = (hConsoleWnd == GetForegroundWindow());

    /*printf("hConsoleWnd = %ld, hDesktopWnd = %ld\n", (long)hConsoleWnd, (long)hDesktopWnd );  */

    if ( !hConsoleWnd || !hDesktopWnd )
        return (FALSE); /* failed */
    if ( !IsWindowVisible(hConsoleWnd) )
        return (FALSE); /* failed */


    /* we will create graphics window of the same size and position as our console window */
    /* to do that we need to get console app window size and position. */
    ret  = GetWindowRect( hConsoleWnd, &rc );
    /*printf( "ConsoleWnd: ret=%d, rc=%ld %ld %ld %ld\n", ret, rc.left, rc.top, rc.right, rc.bottom); */
    /* full screen: "ConsoleWnd: ret=1, rc=-32000 -32000 -31840 -31976" */
    
    ret2 = GetWindowRect( hDesktopWnd, &rc2 );
    /*printf( "DesktopWnd: ret=%d, rc2=%ld %ld %ld %ld\n", ret, rc2.left, rc2.top, rc2.right, rc2.bottom); */
    /* full screen: "DesktopWnd: ret=1, rc2=0 0 1280 1024" */
    
    if ( !(hDesktopDC = GetWindowDC(hDesktopWnd) ) )
        return (FALSE);
    bRectVisible = RectVisible( hDesktopDC, &rc );
    /*printf( "Console rect visible=%d\n", bRectVisible); */
    /* full screen: "Console rect visible=1" */
    
    /*ret3 = GetClipBox(hDesktopDC, &rc3); */
    /*printf( "DesktopClip: ret=%d, rc3=%ld %ld %ld %ld\n", ret3, rc3.left, rc3.top, rc3.right, rc3.bottom); */
    /* full screen: "DesktopClip: ret=3, rc3=0 0 0 0", 3=COMPLEXREGION */
    
    ret3 = ReleaseDC( hDesktopWnd, hDesktopDC );
    if ( !bRectVisible )
        return (FALSE); /* usually happens in MS-DOS full-screen mode, but the API call may fail and return TRUE */

    if ( !IntersectRect(&rc3, &rc2, &rc) ) {
        /*printf( "Console rect invisible\n", bRectVisible); */
        /* full screen: "Console rect invisible" */
        return (FALSE); /* usually happens in MS-DOS full-screen mode */
    }
    
    hInstance = GetModuleHandle (NULL); /* or =(HINSTANCE)GetWindowLong(hConsoleWnd, GWL_HINSTANCE); */
    if ( !dp->pdp->rcPict[2] || !dp->pdp->rcPict[3] ) {
        dp->pdp->rcPict[0] = rc.left;
        dp->pdp->rcPict[1] = rc.top;
        dp->pdp->rcPict[2] = rc.right-rc.left;
        dp->pdp->rcPict[3] = rc.bottom-rc.top;
    }

    /* Fill in window class structure with parameters that describe */
    /* the main window. */
    wc.style         = CS_HREDRAW | CS_VREDRAW;
    wc.lpfnWndProc   = (WNDPROC)WndProcDisplayInputStructure; /* window procedure */
    wc.cbClsExtra    = 0;
    wc.cbWndExtra    = 0;
    wc.hInstance     = hInstance;
    wc.hIcon         = LoadIcon (NULL, IDI_APPLICATION);
    wc.hCursor       = LoadCursor(NULL, IDC_ARROW);
    wc.hbrBackground = (HBRUSH)GetStockObject(WHITE_BRUSH); /*(COLOR_WINDOW+1); */
    wc.lpszMenuName  = NULL;
    wc.lpszClassName = szWindowClassName;
    /* Register the window class and return success/failure code. */
    if ( !RegisterClass(&wc) )
        return FALSE;

    hWnd = CreateWindow(szWindowClassName, szOutputString, WS_OVERLAPPEDWINDOW,
        dp->pdp->rcPict[0], dp->pdp->rcPict[1], dp->pdp->rcPict[2], dp->pdp->rcPict[3],
        /*rc.left, rc.top, rc.right-rc.left, rc.bottom-rc.top, */
        hConsoleWnd, NULL, hInstance, NULL);

    if (!hWnd) {
        /* avoid nasty messages about exception 0x5 in Windows dlls. */
        UnregisterClass(szWindowClassName,    /* address of class name string */
                        hInstance             /* handle of application instance */
                       );
        return (FALSE);
    }

    /* provide window procedure with the pointer to the chemical structure */
    memset( &WinData, 0, sizeof(WinData) );
    WinData.at0   =(inp_ATOM *)inchi_calloc(num_at+1, sizeof(WinData.at0[0]));
    WinData.at1   =(inp_ATOM *)inchi_calloc(num_at+1, sizeof(WinData.at1[0]));

    if ( dp->nEquLabels && dp->nNumEquSets ) {
        WinData.nEquLabels   = (AT_NUMB *)inchi_calloc( num_at+1, sizeof(WinData.nEquLabels[0]));
        WinData.nNumEquSets  = dp->nNumEquSets;
        WinData.nCurEquLabel = 0;
        WinData.nNewEquLabel = dp->nCurEquLabel;
    }

    WinData.nFontSize = dp->sdp.nFontSize;
    WinData.szTitle   = NULL;
    /*WinData.szTitle   = _strdup(szOutputString);*/ /* for testing INCHI_LIB printing */
    if ( inf_at ) {
        DuplicateInfoAtomData( &WinData.inf_at_data, inf_at_data);
    }
    if ( !WinData.at0 || !WinData.at1 || inf_at && !WinData.inf_at_data.at
        || dp->nEquLabels && dp->nNumEquSets && !WinData.nEquLabels
        ) {
        FreeWinData( &WinData );
    } else {
        memcpy( WinData.at0, at, sizeof(at[0])*num_at );
        if ( inf_at )
            memcpy( WinData.inf_at_data.at, inf_at, sizeof(inf_at[0])*num_at );
        if ( WinData.nEquLabels ) {
            memcpy( WinData.nEquLabels, dp->nEquLabels, num_at*sizeof(WinData.nEquLabels[0]));
        }

        memcpy( WinData.at1,  WinData.at0, sizeof(at[0])*num_at );

        WinData.num_at      = num_at;
        WinData.bOrigAtom   = dp->sdp.bOrigAtom;
        WinData.nTimerId    = 0;
        WinData.ulDisplTime = dp->sdp.ulDisplTime;
        if ( dp->sdp.tdp ) {
            WinData.tdp = *dp->sdp.tdp;
        }
    }

    SetWindowLong( hWnd, GWL_USERDATA, (long)&WinData );

    ShowWindow(hWnd, SW_SHOWNORMAL /*SW_SHOW*/);
    UpdateWindow(hWnd);
    
    
    /* Message Loop for Display Canon Struct Window */
    while( IsWindow(hWnd) ) {
        if ( PeekMessage( 
            &msg,            /* pointer to structure for message */
            hWnd,           /* or NULL,*/ /* handle to window */
            0,              /* UINT wMsgFilterMin, */ /* first message */
            0,              /* UINT wMsgFilterMax, */ /* last message */
            PM_REMOVE       /* UINT wRemoveMsg     */ /* removal flags */
           ) ) {
            TranslateMessage( &msg );
            DispatchMessage( &msg );
        } else {
            if ( bSetForeground ) {
                ReallySetForegroundWindow( hWnd );
                bSetForeground = 0;
            }
            SleepEx( 10L, TRUE ); /* provides a nice behavior of the app */
        }
    }
    /* deallocate memory */
    FreeWinData( &WinData );
    /* show console window on the top upon closing hWnd */
    ReallySetForegroundWindow( hConsoleWnd );
    /* avoid nasty messages about exception 0x5 in Windows dlls. */
    UnregisterClass(szWindowClassName,    /* address of class name string */
                    hInstance             /* handle of application instance */
                   );
    /* Save window size and position */
    if ( WinData.rc.bottom > WinData.rc.top && WinData.rc.right > WinData.rc.left ) {
        dp->pdp->rcPict[0] = WinData.rc.left;
        dp->pdp->rcPict[1] = WinData.rc.top;
        dp->pdp->rcPict[2] = WinData.rc.right-WinData.rc.left;
        dp->pdp->rcPict[3] = WinData.rc.bottom-WinData.rc.top;
    }
    if ( WinData.bEsc ) {
        dp->rdp.bEsc = 1;
    }
    return WinData.bEsc? 27:1;

}
#endif
/****************************************************************************/
void MySleep( unsigned long ms )
{
    Sleep( ms );
}

#endif /* } _WIN32 */
#else
int dummyDispStru_c; /* make translation unit non-empty for ANSI-C compatibility */
#endif /* } INCHI_ANSI_ONLY */
