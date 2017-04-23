/*
For now most of the code comes from http://home.hiwaay.net/~taylorc/toolbox/geography/geoutm.html till i actually bother to understand the maths behind it. but finding instructions is hard
*/


var notImpl = function(method="unknown Method"){ var method_ = method; return function(){console.log("not implemented: "+method_); }}

var UTMOptions = { //initial values for WGS84
    equatorialRadius:6378137.0,
    polarRadius: 6356752.3142,
    flattening: 1/298.257223563
};

//public methods will be assigned in a bit
var LatLonToUTM = notImpl("LatLonToUTM");
var UTMToLatLon = notImpl("UTMToLatLon");

(function(){
    var UTMScaleFactor = 0.9996;
    
    //zone - UTM zone to be used for calculating values for x and y. If zone is less than 1 or greater than 60, the routine will determine the appropriate zone from the value of lon.
    //returns [x,y,zone,isSouthernHemisphere]
    LatLonToUTM = function(lat, lon, targetZone = 0)
    {
        var xy = [0,0];
        
        if(targetZone < 1 || targetZone >60){
            targetZone = UTMZoneFromLatLon(lat,lon);
        }
        
        var lambda0 = UTMCentralMeridian (targetZone);
        MapLatLonToXY (DegToRad(lat), DegToRad(lon), lambda0, xy);

        /* Adjust easting and northing for UTM system. */
        xy[0] = xy[0] * UTMScaleFactor + 500000.0;
        xy[1] = xy[1] * UTMScaleFactor;
        if (xy[1] < 0.0){
            xy[1] = xy[1] + 10000000.0;
        }

        retval = [xy[0],xy[1],targetZone, lat < 0];
        return retval;
    }
    
    UTMToLatLon = function(x, y, zone, southernHemisphere)
    {
        var cmeridian;
        	
        x -= 500000.0;
        x /= UTMScaleFactor;
        	
        /* If in southern hemisphere, adjust y accordingly. */
        if (southernHemisphere){
            y -= 10000000.0;
        }
        y /= UTMScaleFactor;
        
        cmeridian = UTMCentralMeridian (zone);
        MapXYToLatLon (x, y, cmeridian);
        	
        return;
    }
    
    function UTMZoneFromLatLon(lat,lon){
        utmZone = 1 + Math.floor((lon + 180) / 6);

        if(lat >= 56.0 && lat < 64.0 && lon >= 3.0 && lon < 12.0){
            utmZone = 32;
        }
        // Special zones for Svalbard
        if(lat >= 72.0 && lat < 84.0) {
            if(lon >= 0.0 && lon < 9.0) {
                utmZone = 31;
            }else if(lon >= 9.0 && lon < 21.0){
                utmZone = 33;
            }else if(lon >= 21.0 && lon < 33.0){
                utmZone = 35;
            }else if(lon >= 33.0 && lon < 42.0) {
                utmZone = 37;
            }
        }
        
        return utmZone;
    }
    
    /*
    * UTMCentralMeridian
    *
    * Determines the central meridian for the given UTM zone.
    *
    * Inputs:
    *     zone - An integer value designating the UTM zone, range [1,60].
    *
    * Returns:
    *   The central meridian for the given UTM zone, in radians, or zero
    *   if the UTM zone parameter is outside the range [1,60].
    *   Range of the central meridian is the radian equivalent of [-177,+177].
    *
    */
    function UTMCentralMeridian (zone)
    {
        var cmeridian;
        cmeridian = DegToRad (-183.0 + (zone * 6.0));
        return cmeridian;
    }
    
    
    /*
    * MapLatLonToXY
    *
    * Converts a latitude/longitude pair to x and y coordinates in the
    * Transverse Mercator projection.  Note that Transverse Mercator is not
    * the same as UTM; a scale factor is required to convert between them.
    *
    * Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
    * GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994.
    *
    * Inputs:
    *    phi - Latitude of the point, in radians.
    *    lambda - Longitude of the point, in radians.
    *    lambda0 - Longitude of the central meridian to be used, in radians.
    *
    * Outputs:
    *    xy - A 2-element array containing the x and y coordinates
    *         of the computed point.
    *
    * Returns:
    *    The function does not return a value.
    *
    */
    function MapLatLonToXY (phi, lambda, lambda0, xy)
    {
        var N, nu2, ep2, t, t2, l;
        var l3coef, l4coef, l5coef, l6coef, l7coef, l8coef;
        var tmp;

        /* Precalculate ep2 */
        ep2 = (Math.pow (UTMOptions.equatorialRadius, 2.0) - Math.pow (UTMOptions.polarRadius, 2.0)) / Math.pow (UTMOptions.polarRadius, 2.0);
    
        /* Precalculate nu2 */
        nu2 = ep2 * Math.pow (Math.cos (phi), 2.0);
    
        /* Precalculate N */
        N = Math.pow (UTMOptions.equatorialRadius, 2.0) / (UTMOptions.polarRadius * Math.sqrt (1 + nu2));
    
        /* Precalculate t */
        t = Math.tan (phi);
        t2 = t * t;
        tmp = (t2 * t2 * t2) - Math.pow (t, 6.0);

        /* Precalculate l */
        l = lambda - lambda0;
    
        /* Precalculate coefficients for l**n in the equations below
           so a normal human being can read the expressions for easting
           and northing
           -- l**1 and l**2 have coefficients of 1.0 */
        l3coef = 1.0 - t2 + nu2;
    
        l4coef = 5.0 - t2 + 9 * nu2 + 4.0 * (nu2 * nu2);
    
        l5coef = 5.0 - 18.0 * t2 + (t2 * t2) + 14.0 * nu2
            - 58.0 * t2 * nu2;
    
        l6coef = 61.0 - 58.0 * t2 + (t2 * t2) + 270.0 * nu2
            - 330.0 * t2 * nu2;
    
        l7coef = 61.0 - 479.0 * t2 + 179.0 * (t2 * t2) - (t2 * t2 * t2);
    
        l8coef = 1385.0 - 3111.0 * t2 + 543.0 * (t2 * t2) - (t2 * t2 * t2);
    
        /* Calculate easting (x) */
        xy[0] = N * Math.cos (phi) * l
            + (N / 6.0 * Math.pow (Math.cos (phi), 3.0) * l3coef * Math.pow (l, 3.0))
            + (N / 120.0 * Math.pow (Math.cos (phi), 5.0) * l5coef * Math.pow (l, 5.0))
            + (N / 5040.0 * Math.pow (Math.cos (phi), 7.0) * l7coef * Math.pow (l, 7.0));
    
        /* Calculate northing (y) */
        xy[1] = ArcLengthOfMeridian (phi)
            + (t / 2.0 * N * Math.pow (Math.cos (phi), 2.0) * Math.pow (l, 2.0))
            + (t / 24.0 * N * Math.pow (Math.cos (phi), 4.0) * l4coef * Math.pow (l, 4.0))
            + (t / 720.0 * N * Math.pow (Math.cos (phi), 6.0) * l6coef * Math.pow (l, 6.0))
            + (t / 40320.0 * N * Math.pow (Math.cos (phi), 8.0) * l8coef * Math.pow (l, 8.0));
    
        return;
    }
    
    /*
    * MapXYToLatLon
    *
    * Converts x and y coordinates in the Transverse Mercator projection to
    * a latitude/longitude pair.  Note that Transverse Mercator is not
    * the same as UTM; a scale factor is required to convert between them.
    *
    * Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
    *   GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994.
    *
    * Inputs:
    *   x - The easting of the point, in meters.
    *   y - The northing of the point, in meters.
    *   lambda0 - Longitude of the central meridian to be used, in radians.
    *
    *
    * Returns:
    *   [lat,lon]
    *
    * Remarks:
    *   The local variables Nf, nuf2, tf, and tf2 serve the same purpose as
    *   N, nu2, t, and t2 in MapLatLonToXY, but they are computed with respect
    *   to the footpoint latitude phif.
    *
    *   x1frac, x2frac, x2poly, x3poly, etc. are to enhance readability and
    *   to optimize computations.
    *
    */
    function MapXYToLatLon (x, y, lambda0)
    {
        /* Get the value of phif, the footpoint latitude. */
        var phif = FootpointLatitude (y);
        	
        /* Precalculate ep2 */
        var ep2 = (Math.pow (UTMOptions.equatorialRadius, 2.0) - Math.pow (UTMOptions.polarRadius, 2.0))
              / Math.pow (UTMOptions.polarRadius, 2.0);
        	
        /* Precalculate cos (phif) */
        var cf = Math.cos (phif);
        	
        /* Precalculate nuf2 */
        var nuf2 = ep2 * Math.pow (cf, 2.0);
        	
        /* Precalculate Nf and initialize Nfpow */
        var Nf = Math.pow (UTMOptions.equatorialRadius, 2.0) / (UTMOptions.polarRadius * Math.sqrt (1 + nuf2));
        var Nfpow = Nf;
        	
        /* Precalculate tf */
        var tf = Math.tan (phif);
        var tf2 = tf * tf;
        var tf4 = tf2 * tf2;
        
        /* Precalculate fractional coefficients for x**n in the equations
           below to simplify the expressions for latitude and longitude. */
        var xnFrac = [];
        var nomin = [1.0,tf];
        var mult2 = [cf,1.0];
        var i,mult1,Nfpow1;
        for(i =0,mult1 = 1,Nfpow1 = Nf; i< 8;i++){
            xnFrac[i] = nomin[i%2]/(mult1*Nfpow1 *mult2[i%2]);
            mult1 *= i+2;
            Nfpow1 *= Nf;
        }
           
        var x1frac = 1.0 / (Nfpow * cf);
        
        Nfpow *= Nf;   /* now equals Nf**2) */
        var x2frac = tf / (2.0 * Nfpow);
        
        Nfpow *= Nf;   /* now equals Nf**3) */
        var x3frac = 1.0 / (6.0 * Nfpow * cf);
        
        Nfpow *= Nf;   /* now equals Nf**4) */
        var x4frac = tf / (24.0 * Nfpow);
        
        Nfpow *= Nf;   /* now equals Nf**5) */
        var x5frac = 1.0 / (120.0 * Nfpow * cf);
        
        Nfpow *= Nf;   /* now equals Nf**6) */
        var x6frac = tf / (720.0 * Nfpow);
        
        Nfpow *= Nf;   /* now equals Nf**7) */
        var x7frac = 1.0 / (5040.0 * Nfpow * cf);
        
        Nfpow *= Nf;   /* now equals Nf**8) */
        var x8frac = tf / (40320.0 * Nfpow);
        
        console.log("same: "+(xnFrac[0] == x1frac));
        console.log("same: "+(xnFrac[1] == x2frac));
        console.log("same: "+(xnFrac[2] == x3frac));
        console.log("same: "+(xnFrac[3] == x4frac));
        console.log("same: "+(xnFrac[4] == x5frac));
        console.log("same: "+(xnFrac[5] == x6frac));
        console.log("same: "+(xnFrac[6] == x7frac));
        console.log("same: "+(xnFrac[7] == x8frac));
        
        /* Precalculate polynomial coefficients for x**n.
           -- x**1 does not have a polynomial coefficient. */
        var x2poly = -1.0 - nuf2;
        var x3poly = -1.0 - 2 * tf2 - nuf2;
        var x4poly = 5.0 + 3.0 * tf2 + 6.0 * nuf2 - 6.0 * tf2 * nuf2 - 3.0 * (nuf2 * nuf2) - 9.0 * tf2 * (nuf2 * nuf2);
        var x5poly = 5.0 + 28.0 * tf2 + 24.0 * tf4 + 6.0 * nuf2 + 8.0 * tf2 * nuf2;
        var x6poly = -61.0 - 90.0 * tf2 - 45.0 * tf4 - 107.0 * nuf2 + 162.0 * tf2 * nuf2;
        var x7poly = -61.0 - 662.0 * tf2 - 1320.0 * tf4 - 720.0 * (tf4 * tf2);
        var x8poly = 1385.0 + 3633.0 * tf2 + 4095.0 * tf4 + 1575 * (tf4 * tf2);
        	
        /* Calculate latitude */
        var lat = phif + x2frac * x2poly * (x * x)
        	+ x4frac * x4poly * Math.pow (x, 4.0)
        	+ x6frac * x6poly * Math.pow (x, 6.0)
        	+ x8frac * x8poly * Math.pow (x, 8.0);
        	
        /* Calculate longitude */
        var lon = lambda0 + x1frac * x
        	+ x3frac * x3poly * Math.pow (x, 3.0)
        	+ x5frac * x5poly * Math.pow (x, 5.0)
        	+ x7frac * x7poly * Math.pow (x, 7.0);
        	
        return [lat,lon];
    }

    
    /*
    * DegToRad
    *
    * Converts degrees to radians.
    *
    */
    function DegToRad (deg)
    {
        return (deg / 180.0 * Math.PI)
    }
    
    /*
    * RadToDeg
    *
    * Converts radians to degrees.
    *
    */
    function RadToDeg (rad)
    {
        return (rad / Math.PI * 180.0)
    }

    /*
    * ArcLengthOfMeridian
    *
    * Computes the ellipsoidal distance from the equator to a point at a
    * given latitude.
    *
    * Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
    * GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994.
    *
    * Inputs:
    *     phi - Latitude of the point, in radians.
    *
    * Globals:
    *     UTMOptions.equatorialRadius - Ellipsoid model major axis.
    *     UTMOptions.polarRadius - Ellipsoid model minor axis.
    *
    * Returns:
    *     The ellipsoidal distance of the point from the equator, in meters.
    *
    */
    function ArcLengthOfMeridian (phi)
    {
        /* Precalculate n */
        var n = (UTMOptions.equatorialRadius - UTMOptions.polarRadius) / (UTMOptions.equatorialRadius + UTMOptions.polarRadius);

        /* Precalculate alpha */
        var alpha = ((UTMOptions.equatorialRadius + UTMOptions.polarRadius) / 2.0)
           * (1.0 + (Math.pow (n, 2.0) / 4.0) + (Math.pow (n, 4.0) / 64.0));

        /* Precalculate beta */
        var beta = (-3.0 * n / 2.0) + (9.0 * Math.pow (n, 3.0) / 16.0)
           + (-3.0 * Math.pow (n, 5.0) / 32.0);

        /* Precalculate gamma */
        var gamma = (15.0 * Math.pow (n, 2.0) / 16.0)
            + (-15.0 * Math.pow (n, 4.0) / 32.0);
    
        /* Precalculate delta */
        var delta = (-35.0 * Math.pow (n, 3.0) / 48.0)
            + (105.0 * Math.pow (n, 5.0) / 256.0);
    
        /* Precalculate epsilon */
        var epsilon = (315.0 * Math.pow (n, 4.0) / 512.0);
    
    /* Now calculate the sum of the series and return */
        var result = alpha
        * (phi + (beta * Math.sin (2.0 * phi))
            + (gamma * Math.sin (4.0 * phi))
            + (delta * Math.sin (6.0 * phi))
            + (epsilon * Math.sin (8.0 * phi)));

    return result;
    }
    
    /*
    * FootpointLatitude
    *
    * Computes the footpoint latitude for use in converting transverse
    * Mercator coordinates to ellipsoidal coordinates.
    *
    * Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
    *   GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994.
    *
    * Inputs:
    *   y - The UTM northing coordinate, in meters.
    *
    * Returns:
    *   The footpoint latitude, in radians.
    *
    */
    function FootpointLatitude (y)
    {        
        /* Precalculate n (Eq. 10.18) */
        var n = (UTMOptions.equatorialRadius - UTMOptions.polarRadius) / (UTMOptions.equatorialRadius + UTMOptions.polarRadius);
        	
        /* Precalculate alpha_ (Eq. 10.22) */
        /* (Same as alpha in Eq. 10.17) */
        var alpha_ = ((UTMOptions.equatorialRadius + UTMOptions.polarRadius) / 2.0)
            * (1 + (Math.pow (n, 2.0) / 4) + (Math.pow (n, 4.0) / 64));
        
        /* Precalculate y_ (Eq. 10.23) */
        var y_ = y / alpha_;
        
        /* Precalculate beta_ (Eq. 10.22) */
        var beta_ = (3.0 * n / 2.0) + (-27.0 * Math.pow (n, 3.0) / 32.0)
            + (269.0 * Math.pow (n, 5.0) / 512.0);
        
        /* Precalculate gamma_ (Eq. 10.22) */
        var gamma_ = (21.0 * Math.pow (n, 2.0) / 16.0)
            + (-55.0 * Math.pow (n, 4.0) / 32.0);
        	
        /* Precalculate delta_ (Eq. 10.22) */
        var delta_ = (151.0 * Math.pow (n, 3.0) / 96.0)
            + (-417.0 * Math.pow (n, 5.0) / 128.0);
        	
        /* Precalculate epsilon_ (Eq. 10.22) */
        var epsilon_ = (1097.0 * Math.pow (n, 4.0) / 512.0);
        	
        /* Now calculate the sum of the series (Eq. 10.21) */
        var result = y_ + (beta_ * Math.sin (2.0 * y_))
            + (gamma_ * Math.sin (4.0 * y_))
            + (delta_ * Math.sin (6.0 * y_))
            + (epsilon_ * Math.sin (8.0 * y_));
        
        return result;
    }
})();
