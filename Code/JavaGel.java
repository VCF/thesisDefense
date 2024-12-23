/*
  JavaGel
  Charles Tilford		11/11/97
*/

import java.awt.*;
import java.awt.image.*;
import java.applet.*;
import java.util.*;
import java.net.*;
import java.io.*;

class Band {
    int x;
    int y;
    int column;
    int plate;
    Band (int a,int b,int c, int p) {
        x = a;
        y = b;
        column = c;
        plate = p;
    }
}
	
public class JavaGel extends Applet {
    int xloc = 0;		// X position of gph
    int yloc = 32;		// Y position of gph
    int xwid = 5;		// X width of band box
    int ywid = 5;		// Y width of band box
    int MouseX = 0;		// X position of mouse
    int MouseY = 0;		// Y position of mouse
    int gphWidth, gphHeight;
    int width = 520, height = 1000;
    int xinit = -1, yinit = -1, xbox, ybox;
    int plateNum, $wells, tiers;
    Image gph, offscreenImage;
    Graphics osg;
    String msg = "", debug = "";
    Color red = new Color(255, 0, 0);
    Color green = new Color(0, 255, 0);
    Color blue = new Color(0, 0, 255);
    Color yellow = new Color(255, 255, 0);
    Color black = new Color(0, 0, 0);
    Color white = new Color(255, 255, 255);
    Vector data = new Vector(10,2);
    Vector undoData = new Vector(10,2);
    char sep = 's';
    char sep2 = 'q';
    int moveWell = -1, newband = 0;
    int wells[][];
    int qual[] = {-1,-1};		// Holds Quality values for Top and Bottom frame
    int temp[] = {-1,-1};		// Same for Temp.
    int lowT = 54, topT = 64;	// Range of temps to display on screen
    int targ = 11;				// Size of target box for Qual / Temp
    int targsep = 2;			// Separation between target boxes
    String qualities[]	= { "Excellent", "Reasonable", "Poor", "Useless" };
    String qualCodes[]	= { "A",         "B",          "C",    "D" };
    int navx = 3*width/5;		// navx defines left boundary of clickable area.
    int navx2 = navx + qualities.length*(targ+targsep)+5;		
    int navy[] = {3,targ+3+3};	// navy defines top boundaries
	
    String workKeys[] = { "W", "S", "A", "D", "Z", "PU"	 };
    int keyNums[]     = { 119, 115, 97,  100, 122, Event.PGUP };
    String workModes[]	= { "Add Wells", "Search", "Add Bands", "Delete", "Undo Last", "Submit"	};
    int dMode = 0, oMode = 0, uMode = 0;
    boolean undoActive = true, mouseClicked = false, noSubmit = false, hideBox = false;
    boolean controls[] = new boolean[97];
    protected Frame frame;
    boolean compareTier[][] = new boolean[2][97];
    String rout[] = {	"mouseDown",	"mouseDrag",	"mouseUp",
                        "findBands",	"intensity",	"findAllColumns",
                        "findColumn",	"allBand",		"keyDown",
                        "paint", "Total Time" };
    long tim[] = new long[ rout.length +1 ];
    long subTim;
    boolean deBug = false;
	
    // + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 	
    public void init() {
        gph = getImage( getDocumentBase(), getParameter("fname") );
        offscreenImage = createImage(size().width, size().height);
        osg = offscreenImage.getGraphics();
        tiers = Integer.parseInt( getParameter("tiers") );
        $wells = (tiers+1) * 2;
        wells = new int[$wells][2];
        for (int i=0; i < $wells; i++) {
            wells[i][0] = 0;
            wells[i][1] = 0;
        }
        controls[1] = true; controls[6] = true; controls[10] = true; controls[96] = true; 
        controls[2] = true; controls[95] = true;
        String parseBands = getParameter("bandstuff");
        String parseWells = getParameter("wellstuff");
        String hybStats = getParameter("hybstats");
        int j=-1, x=0, y=0, w=0; String tmp;
		
        if (parseWells.length() > 0) {
            for (int i=0; i < parseWells.length(); i++) {
                if (w >= $wells) break;
                if (parseWells.charAt(i) == sep) {
                    tmp = parseWells.substring(j+1, i);	// Finds the stuff between the separator.
                    j = i;
                    if (x == 0) {
                        x = Integer.parseInt(tmp); 
                    } else {
                        y = Integer.parseInt(tmp);
                        if (x > size().width) x = size().width;
                        if (y > size().height - yloc) y = size().height - yloc;
                        wells[w][0] = x;
                        wells[w][1] = y;
                        x = y = 0; w++;
						
                    }
                }
            }
        }
		
        j = -1;
        if (parseBands.length() > 0) {
            for (int i=0; i < parseBands.length(); i++) {
                if (parseBands.charAt(i) == sep) {
                    tmp = parseBands.substring(j+1, i);	// Finds the stuff between the separator.
                    j = i;
                    if (x == 0) {
                        x = Integer.parseInt(tmp); 
                    } else {
                        y = Integer.parseInt(tmp);
                        if (x > size().width) x = size().width;
                        if (y > size().height - yloc) y = size().height - yloc;
                        addBand (x,y);
                        x = y = 0;
                    }
                }
            }
        }

        j = -1; w=0;
        if (hybStats.length() > 0) {
            for (int i=0; i < hybStats.length(); i++) {
                if (hybStats.charAt(i) == sep2) {
                    tmp = hybStats.substring(j+1, i);	// Finds the stuff between the separator.
                    j = i;
                    if (x == 0) {
                        for (int k = 0; k < qualities.length; k++) {
                            if (tmp.equalsIgnoreCase(qualities[k])) { qual[w] = k; }
                        }
                        x = 1; 
                    } else {
                        temp[w] = Integer.parseInt(tmp);
                        x = 0; w++;
						
                    }
                }
            }
        }
		
        // Taken from http://www.vivids.com/java/AllHtmlSource/MultiText.java.html...
        Component c = this;
        frame = null;
        while(c != null) {
            c = c.getParent();
            if(c instanceof Frame) {
                frame = (Frame)c;
                break;
            }
        }
        //frame.setCursor(Frame.CROSSHAIR_CURSOR);
        // ...
        findAllColumns ();
    }
    // + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 	
    public void update(Graphics g) {
        paint(g);
    }	
    //	Mouse commands - page 523
    // + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 	
    public boolean mouseDown(Event evtObj, int x, int y) {
        Date t1 = new Date(); msg = "";
        hideBox = false;
        OUT: {
            if (y < yloc) {				// Click in menu region
                mouseClicked = false; xbox = -1; ybox = -1;
                int xd, i = -1;
                if ( (y >= navy[0]) && (y <= navy[0]+targ) ) i = 0;  // Check vertical position of click
                if ( (y >= navy[1]) && (y <= navy[1]+targ) && (tiers == 3) ) i = 1;
                if (i == -1) break OUT;
                if (x >= navx) {		// Click to right of quality boundary
                    if (x >= navx2) {	// Click to right of temperature boundary
                        for (int k = lowT; k <= topT; k++) {
                            xd = navx2 + (k-lowT)*(targ + targsep);
                            if ( (x >= xd) && (x <= xd+targ) ) {
                                temp[i] = k;
                                break OUT;
                            }
                        }
                        break OUT;
                    }
                    for (int k = 0; k < qualCodes.length; k++) {	// Draw in quality codes
                        xd = navx + k*(targ + targsep);
                        if ( (x >= xd) && (x <= xd+targ) ) {
                            qual[i] = k;
                        }
                    }
                }
                break OUT;
            }
            mouseClicked = true;
            xinit = x - xloc; yinit = y - yloc;
            for (int i=0; i < $wells; i++) {
                if ( (xinit >= wells[i][0]-xwid) && (xinit <= wells[i][0]+xwid) && (yinit >= wells[i][1]-ywid) && (yinit <= wells[i][1]+ywid) ) {
                    moveWell = i;
                    break;
                }
            }
        }
        Date t2 = new Date();
        tim[0] += t2.getTime() - t1.getTime();
        repaint ();
        return true;
    }
    // + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 	
    public boolean mouseDrag(Event evtObj, int x, int y) {
        Date t1 = new Date();
        long otherTim = 0;
        OUT: {
            if (moveWell >= 0) {
                wells[moveWell][0] = x - xloc;
                wells[moveWell][1] = y - yloc;
                findAllColumns ();
                otherTim += subTim;
                break OUT;
            }
            if (dMode == 0) return true;
            xbox = x - xloc; ybox = y - yloc;
        }
        Date t2 = new Date();
        tim[1] += t2.getTime() - t1.getTime() - otherTim;
        repaint ();
        return true;
    }
    // + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 	
    public boolean mouseUp(Event evtObj, int x, int y) {
        Date t1 = new Date(); msg = "";
        long otherTim = 0;
        boolean drawNow = false;
        OUT: {
            if (!mouseClicked) break OUT;
            mouseClicked = false;
            x -= xloc; y -= yloc;
            if (moveWell >= 0) {
                wells[moveWell][0] = x;
                wells[moveWell][1] = y;
                findAllColumns ();
                otherTim += subTim;
                moveWell = -1; drawNow = true;
                break OUT;
            }
            Band da = new Band(-1,-1,-1,-1);
            if (xbox < 0) { xbox = xinit; ybox = yinit; }
            boolean moreWork = true;
            if ((dMode == 3) || ( (dMode !=0) && (xbox != xinit) ) ) {		// Band Removal
                undoData.removeAllElements(); uMode = 3;
                msg = "Bands Removed: ";
                int a = xinit; int b = yinit; int c = xbox; int d = ybox;
                if (a > xbox) { a = xbox; c = xinit; }
                if (b > ybox) { b = ybox; d = yinit; }
                for (int i=0; i < data.size(); i++) {
                    da = (Band) data.elementAt(i);
                    if ( (a < da.x+xwid) && (c > da.x-xwid) && (b < da.y+ywid) && (d > da.y-ywid) ) {
                        undoData.addElement( (Band) data.elementAt(i) );
                        msg += da.column + ", ";
                        compareTier[da.plate][da.column] = false;
                        data.removeElementAt(i);
                        i--;
                    }
                }
                xbox = -1; ybox = -1; drawNow = true;
                break OUT;
            }
            if (dMode == 0) {
                for (int i=0; i < $wells; i++) {
                    if ( (wells[i][0] == 0) && (wells[i][1] == 0) ) {
                        wells[i][0] = x;
                        wells[i][1] = y;
                        moreWork = false;
                        break;
                    }
                }
                if (moreWork) {  // Find the closest well to assign here
                    float min = 100000, d = 0, dx = 0, dy = 0;
                    int j = 0;
                    for (int i=0; i < $wells; i++) {
                        dx = (wells[i][0] - x); dx *= dx;
                        dy = (wells[i][1] - y); dy *= dy;
                        d = (float) Math.sqrt(dx + dy);
                        if (d < min) {
                            min = d; j = i;
                        }
                    }
                    wells[j][0] = x;
                    wells[j][1] = y;
                }
                for (int i=1; i < $wells; i++) {  // Sort by y values
                    int j = i;
                    while ( (j > 0) && (wells[j][1] < wells[j-1][1]) ) {
                        if ( (wells[j][0] == 0) && (wells[j][1] == 0) ) break;
                        x = wells[j][0]; wells[j][0] = wells[j-1][0]; wells[j-1][0] = x;
                        x = wells[j][1]; wells[j][1] = wells[j-1][1]; wells[j-1][1] = x;
                        j--;
                    }
                }
                for (int i=1; i < $wells; i+=2) {  // Sort by x values pairwise
                    if ( wells[i][0] < wells[i-1][0])  {
                        if ( (wells[i][0] == 0) && (wells[i][1] == 0) ) break;
                        x = wells[i][0]; wells[i][0] = wells[i-1][0]; wells[i-1][0] = x;
                        x = wells[i][1]; wells[i][1] = wells[i-1][1]; wells[i-1][1] = x;
                    }
                }
                xwid = (wells[1][0]-wells[0][0])/100;
                if (xwid < 4) xwid = 4;
                ywid = (wells[2][1]-wells[0][1])/14;
                if (ywid < 4) ywid = 4;
				
                findAllColumns ();
                otherTim += subTim;
                msg = "Added well";
            }
            if (dMode == 1) {		// Auto Band Hunt
                undoData.removeAllElements(); uMode = dMode;
                findBands (x, y);
                otherTim += subTim;
            }
            if (dMode == 2) {		// Manual Band Addition
                undoData.removeAllElements(); uMode = dMode;
                addBand (x,y);
                otherTim += subTim;
                msg = "Added new band";
            }
            drawNow = true;
            xbox = -1; ybox = -1;
        }
        Date t2 = new Date();
        tim[2] += t2.getTime() - t1.getTime() - otherTim;
        if (drawNow) repaint();
        return true;
    }
    // + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 	
    void findBands (int a, int b) {
        Date t1 = new Date();
        long otherTim = 0;
        int j = tiers; Date dt = new Date(); long ti = dt.getTime();
        while ( (wells[j*2][1] > b) && (wells[j*2+1][1] > b) && (j>0) || (wells[j*2+1][1] == 0)) j--;
        float m = wells[j*2][1] - wells[j*2+1][1]; float c = a; float d = wells[j*2][1];
        if (m != 0) {
            m = (wells[j*2][0] - wells[j*2+1][0]) / m;
            c = (wells[j*2][0] - (m * wells[j*2][1]) + (m*m*a) + m*b) / (1 + m*m);
            d = m*a - m*c + b;
        }
        int column = (int)Math.floor ( 0.5 + (49 * (c - wells[j*2][0]) / (wells[j*2+1][0] - wells[j*2][0])));
        int x1 = wells[j*2][0]; int x2 = wells[j*2+1][0]; int y1 = wells[j*2][1]; int y2 = wells[j*2+1][1];
        float dx = (float)(x2 - x1) / 49; float dy = (float)(y2 - y1) / 49;
        int Ch = 3; int Cw = (int) (dx / 2) - 1;
        int BGh = Ch; 
        int xr = (int) (x1 + a - c); int yr = (int) (y1 + b - d);	//Reference positions
        int xt = (int) (xr + (dx * column) + 0.5); int yt = (int) (yr + (dy * column) + 0.5);
        int dB = Ch + BGh + 1;
        float band		= intensity (xt - Cw, yt - Ch, Cw * 2 + 1,  Ch * 2 + 1);
        otherTim += subTim;
        float top		= intensity (xt - Cw, yt - dB, Cw * 2 + 1,  BGh);
        otherTim += subTim;
        float bottom	= intensity (xt - Cw, yt + dB, Cw * 2 + 1,  BGh);
        otherTim += subTim;
		
        float StN = band/bottom, signal=0, noise = 0;
        if (top > bottom) { StN = band/top; dB = -dB; }
        int totband = 0; newband = 0;
		
        for (int s=1; s <= 48; s++) {
            int xc = (int) (xr + (dx * s) + 0.5); int yc = (int) (yr + (dy * s) + 0.5);
            signal	= intensity (xc - Cw, yc - Ch, Cw * 2 + 1,  Ch * 2 + 1);
            otherTim += subTim;
            noise	= intensity (xc - Cw, yc + dB, Cw * 2 + 1,  BGh);
            if ( (signal/noise) <= StN ) {
                otherTim += subTim;
                totband++;
                addBand (xc,yc);
                otherTim += subTim;
            }
        }
        ti = 1000 * (ti - dt.getTime()) ;
        msg = "Found a total of " + totband + " bands, " + newband + " being new. (" + (int)ti + " sec).";
        Date t2 = new Date();
        tim[3] += t2.getTime() - t1.getTime() - otherTim;
        subTim = t2.getTime() - t1.getTime();
    }
    // + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 	
    //	a- left		b- top	w- width	h- height
    int intensity (int a, int b, int w, int h) {
        Date t1 = new Date();
        int pixels[] = new int[w*h]; int c = 0;
        PixelGrabber pg = new PixelGrabber(gph, a, b, w, h, pixels, 0, w);
        try {
            pg.grabPixels();
        } catch (InterruptedException e) {
            showStatus("Interrupted waiting for pixels!");
            return -1;
        }
        for (int i=0; i < w * h; i++) {
            c += 0xff & (pixels[i]);
        }
        Date t2 = new Date();
        subTim = t2.getTime() - t1.getTime();
        tim[4] += subTim;
        return c;
    }

    // + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 	
    void findAllColumns () {
        Date t1 = new Date();
        long otherTim = 0;
        Band da = new Band(-1,-1,-1,-1);
        for (int i=0; i < 97; i++) {
            compareTier[0][i] = compareTier[1][i] = false;
        }
		
        for (int i=0; i < data.size(); i++) {
            da = (Band) data.elementAt(i);
            data.setElementAt(new Band (da.x,da.y,findColumn(da.x,da.y), plateNum), i);
            otherTim += subTim;
        }
        Date t2 = new Date();
        tim[5] += t2.getTime() - t1.getTime() - otherTim;
        subTim = t2.getTime() - t1.getTime();
        return;
    }
    // + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 	
    int findColumn (int a, int b) {
        Date t1 = new Date();
        int j = tiers;
        while ( (wells[j*2][1] > b) && (wells[j*2+1][1] > b) && (j>0) || (wells[j*2+1][1] == 0)) j--;
        float m = wells[j*2][1] - wells[j*2+1][1]; float c = a; float d = wells[j*2][1];
        if (m != 0) {
            m = (wells[j*2][0] - wells[j*2+1][0]) / m;
            c = (wells[j*2][0] - (m * wells[j*2][1]) + (m*m*a) + m*b) / (1 + m*m);
            d = m*a - m*c + b;
        }
        int column = (int)Math.floor ( 0.5 + (49 * (c - wells[j*2][0]) / (wells[j*2+1][0] - wells[j*2][0])));
        if (j%2 != 0) column += 48;
        plateNum = (int)Math.floor( j / 2);
        if ((column < 1) || (column > 96)) column = 0;
        compareTier[plateNum][column] = true;
        Date t2 = new Date();
        subTim = t2.getTime() - t1.getTime();
        tim[6] += subTim;
        return column;
    }
    // + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 	
    void addBand (int a, int b) {
        Date t1 = new Date();
        long otherTim = 0;
        int c = findColumn(a,b);
        otherTim += subTim;
        boolean found = false;  Band ada = new Band(-1,-1,c, plateNum);
        for (int i=0; i < data.size(); i++) {
            ada = (Band) data.elementAt(i);
            if ( (a >= ada.x-xwid) && (a <= ada.x+xwid) && (b >= ada.y-ywid) && (b <= ada.y+ywid) ) {
                data.setElementAt(new Band (a,b,c,plateNum), i);
                found = true;
                break;
            }
        }
        if (!found) { data.addElement( new Band (a,b,c,plateNum) ); newband++; }
        if (undoActive) {
            if ( !( found & (dMode == 1) ) ) undoData.addElement( new Band (a,b,c,plateNum) );
        }
        Date t2 = new Date();
        tim[7] += t2.getTime() - t1.getTime() - otherTim;
        subTim = t2.getTime() - t1.getTime();
        return;
    }
    // + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 	
    public boolean keyDown(Event evtObj, int key) {
        Date t1 = new Date();
        long otherTim = 0;
        OUT: {
            if ( (wells[0][1] == 0) || (wells[1][1] == 0) ) break OUT;
            oMode = dMode;
            Band da = new Band(-1,-1,-1,-1);
            Band db = new Band(-1,-1,-1,-1);
            Enumeration en;
			
            if (key == keyNums[0]) {
                dMode = 0;
            } else if (key == keyNums[1]) {
                dMode = 1;
            } else if (key == keyNums[2]) {
                dMode = 2;
            } else if (key == keyNums[3]) {
                dMode = 3;
            } else if (key == 32) {
                hideBox = true;
            } else if (key == 104) {
                deBug = !deBug;
            } else if (key == keyNums[4]) {
                if (undoData.isEmpty()) {
                    showStatus("Undo Buffer is Empty");
                    break OUT;
                }
                en = undoData.elements(); int n = 0;
                undoActive = false;
                if (uMode == 3) {
                    msg = "Re-added ";
                    while (en.hasMoreElements() ) {
                        da = (Band) en.nextElement();
                        showStatus("Re-adding Point at x=" + da.x + " y=" + da.y);
                        addBand (da.x,da.y); n++;
                        otherTim += subTim;
                    }
                    uMode = 1;
                } else {
                    msg = "Removed ";
                    while (en.hasMoreElements() ) {
                        da = (Band) en.nextElement();
                        for (int i=0; i < data.size(); i++) {
                            db = (Band) data.elementAt(i);
                            if ( (da.x == db.x) && (da.y == db.y) ) {
                                data.removeElementAt(i); n++;
                                showStatus("Removing Point at x=" + da.x + " y=" + da.y);
                                i--;
                            }						
                        }
                    }
                    uMode = 3;
                }
                msg += n + " bands from undo buffer.";
                undoActive = true;
                break OUT;
            } else if (key == keyNums[5]) {
                URL theURL = null;
                String url = "http://devo.wi.mit.edu/page-lab/YMAP/RHJavaEnd.cgi?";
                String target = "/page-lab/YMAP/RHJavaEnd.cgi?";
                String wellStr = "wellstuff=", bandStr = "&fyle=" + getParameter("fyle");
                for (int i = 0; i < $wells; i++) {
                    wellStr += Integer.toString(wells[i][0]) + sep + Integer.toString(wells[i][1]) + "'";
                }
                bandStr += "&bandstuff=";
                en = data.elements();
                while (en.hasMoreElements() ) {
                    da = (Band) en.nextElement();
                    bandStr += Integer.toString( da.x ) +		sep + Integer.toString( da.y )		+ sep;
                    bandStr += Integer.toString( da.column )  +	sep + Integer.toString( da.plate )	+ "'";
                }
                target += wellStr + bandStr;
                target += "&quals=";
                for (int i = 0; i <= (tiers-1)/2; i++) {
                    if ( (qual[i] < 0) || (temp[i] < 0) ) {
                        noSubmit = true;
                        msg = "You need to enter quality data before proceeding.";
                        break OUT;
                    }
                    target += qualities[ qual[i] ] + "'" + temp[i] + "'";
                }
                if ($wells < 7) $wells = 5;
                target += "&count=" + $wells;
				
                try { theURL = new URL("http","devo.wi.mit.edu",target); }
                catch ( MalformedURLException e) {
                    showStatus("Bad URL: " + theURL);
                }			
                getAppletContext().showDocument(theURL);
                showStatus(target);
            }
        }
        Date t2 = new Date();
        tim[8] += t2.getTime() - t1.getTime() - otherTim;
        repaint();
        return true;
    }
    // + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 	
    public boolean keyUp(Event evtObj, int key) {
        //		if (dMode == 3) dMode = oMode;
        hideBox = false;
        repaint();
        return true;
    }
    // + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 	
    public void paint(Graphics g) {
        Date t1 = new Date();
        int fType		= 0;
        int fSize[][]	= { {9,10}, {14,16}, {9,11} };		
        osg.setColor(white);
        osg.fillRect(0,0,width,height);
        Font f = osg.getFont();
        f = new Font("TimesRoman",Font.BOLD,9);
        osg.setFont(f);
        FontMetrics fm = osg.getFontMetrics();
        int h = fm.getHeight();
        if (h == 7) fType = 1;	// Check to see if Internet explorer is setting the fonts funny
		
        f = new Font("TimesRoman",Font.BOLD,fSize[0][fType]);
        osg.setFont(f);
        char ch = 'X';
        String tmp = "";
        int w = fm.charWidth(ch);
        int xd = 0, yd = 0, a = 0, b = 0, j = 0, cols = 3, colwidth = 5;
		
        for (int i = 0; i <= (tiers-1)/2; i++) {			// Draw the clickable quality data
            yd = navy[i] + targ-2;
            for (int k = 0; k < qualCodes.length; k++) {	// Draw in quality codes
                xd = navx + k*(targ + targsep);
                if (qual[i] == k) {
                    osg.setColor(yellow);
                    osg.fillRect(xd, navy[i], targ, targ);
                    osg.setColor(black);
                } else {
                    osg.setColor(blue);
                }
                osg.drawString(qualCodes[k], xd+1, yd);
            }

            for (int k = lowT; k <= topT; k++) {			// Draw in temperature options
                xd = navx2 + (k-lowT)*(targ + targsep);
                if (temp[i] == k) {
                    osg.setColor(yellow);
                    osg.fillRect(xd, navy[i], targ, targ);
                    osg.setColor(black);
                    osg.drawString(Integer.toString(k), xd, yd);
                } else if (k%2 == 0) {
                    osg.setColor(green);
                    osg.drawString(Integer.toString(k), xd, yd);
                }
            }
        }
		
        f = new Font("TimesRoman",Font.BOLD,fSize[1][fType]);
        osg.setFont(f);
        h = 14; int compIE = fType * 0;

        for (int i = 0; i < workModes.length; i++) {
            xd = 5 + (i%cols) * width / colwidth;
            yd = (1+(i / cols)) * h;
            osg.setColor(yellow);
            if (dMode == i) osg.fillRect(xd-2, yd - h+2, width / colwidth - 5, h );
            osg.setColor(black);
            osg.drawString(workKeys[i], xd, compIE + yd);
            xd += w*3;
            osg.drawString(workModes[i], xd, compIE + yd);
        }

        osg.drawImage(gph, xloc, yloc, this);
        gphWidth = gph.getWidth(this); gphHeight = gph.getHeight(this); 
        osg.setColor(red);
        for (int i=0; i < $wells; i++) {
            if ( (wells[i][0] != 0) || (wells[i][1] != 0) ) {
                osg.drawRect(wells[i][0]-xwid+xloc, wells[i][1]-ywid+yloc, xwid * 2, ywid * 2);
                xd = i + 1; tmp = Integer.toString(xd);
                osg.drawString(tmp , wells[i][0]-xwid+xloc, wells[i][1]-ywid+yloc-1);
            }
        }
        f = new Font("TimesRoman",Font.BOLD,fSize[2][fType]);
        osg.setFont(f);

        osg.setColor(green);
        if ( (dMode == 0) || (moveWell >= 0) ) {
            for (int i=0; i < $wells; i+=2) {
                if ( (wells[i][1] != 0) && (wells[i+1][1] != 0) ) {
                    float dx = (float)(wells[i+1][0] - wells[i][0]) / 49; float dy = (float)(wells[i+1][1] - wells[i][1]) / 49;
                    for (j=1; j<=48; j++) {
                        xd = ((int)(dx*j) + wells[i][0]+xloc); yd = ((int)(dy*j) + wells[i][1]+yloc);
                        osg.drawLine(xd-1, yd, xd+1, yd);
                        osg.drawLine(xd, yd-1, xd, yd+1);
                    }
                }
            }
        }
		
        if (!hideBox) {
            Enumeration en = data.elements(); Band bd = new Band(-1,-1,-1,-1);
            while (en.hasMoreElements() ) {
                osg.setColor(blue);
                bd = (Band) en.nextElement();
                if ((bd.column > 0) & (bd.column <= 96) & (controls[bd.column])) osg.setColor(green);
                xd = bd.x-xwid+xloc; yd = bd.y-ywid+yloc;
                osg.drawRect(xd, yd, xwid * 2, ywid * 2);
                tmp = Integer.toString(bd.column);
                osg.drawString(tmp , xd, yd-1);
            }
        } else if (tiers == 3) {
            osg.setColor(yellow);
            for (int i=1; i < 97; i++) {
                if (compareTier[0][i] != compareTier[1][i]) {
                    a = 0; b = i;
                    if (i > 48) { a = 2; b -= 48; }
                    for (j=a; j < 7; j+= 4) {
                        xd = (int)(wells[j][0] + b * (wells[j+1][0] - wells[j][0]) / 49) + xloc;
                        yd = (int)(wells[j][1] + b * (wells[j+1][1] - wells[j][1]) / 49) + yloc;
                        int xp[] = {xd, xd-3, xd+3};
                        int yp[] = {yd, yd-3, yd-3};
                        osg.fillPolygon(xp,yp,3);
                    }
                }
            }
        }
		
        if ( (xbox >= 0) && (ybox >= 0) ) {
            osg.setColor(yellow);
            a = xinit; b = yinit; int c = xbox - xinit; int d = ybox - yinit;
            if (a > xbox) { a = xbox; c = -c; }
            if (b > ybox) { b = ybox; d = -d; }
            osg.drawRect(a+xloc, b+yloc, c, d);
        }

        if (noSubmit) {
            f = new Font("TimesRoman",Font.BOLD,36);
            osg.setFont(f);
            osg.setColor(red);
            osg.drawString("Enter Qualities and Temps!", 0, 150);
            noSubmit = false;
        }
        Date t2 = new Date();
        tim[9] += t2.getTime() - t1.getTime();
        if (deBug) qualControl(yloc);
        g.drawImage(offscreenImage, 0, 0, this);
        if (msg.length() > 0) showStatus(msg);
    }
    // + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 	
    void qualControl(int yd) {
        Date t1 = new Date();
        tim[10] = t1.getTime() - tim[10];
        Color orange = new Color(255,128,0);
        osg.setColor(orange);
        for (int k = 0; k < rout.length; k++) {
            yd += 10;
            double per = Math.floor( 100 * tim[k]/tim[10]);
            osg.drawString( rout[k],					50, yd);
            osg.drawString( tim[k] + " msec",			150, yd);
            osg.drawString( Double.toString(per) + "%", 220, yd);
        }
        yd += 12;
        Font f2 = new Font("TimesRoman",Font.BOLD,14);
        osg.setFont(f2);
        osg.drawString( "Hit h to get rid of debugging information.", 50, yd);
        for (int k = 0; k < rout.length; k++) {
            tim[k] = 0;
        }
        tim[10] = t1.getTime(); 
    }
    // + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 	
}
