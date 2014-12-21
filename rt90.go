// Package rt90 transforms swedish grid (RT90) coordinates to WGS84.
// See: http://www.lantmateriet.se/en/Maps-and-geographic-information/GPS-and-geodetic-surveys/Reference-systems/Two-dimensional-systems/RT-90/
package rt90

import "math"

// WGS84FromRT90 transforms RT90 coordinates to WGS84 coordinates.
func WGS84FromRT90(x, y float64) (lat float64, long float64) {
	lat, long = gaussKrüger(x, y)
	return
}

func gaussKrüger(x, y float64) (float64, float64) {
	var axis float64 = 6378137.0
	var flattening float64 = 1.0 / 298.257222101
	var centralMeridian float64 = 15.0 + 48.0/60.0 + 22.624306/3600.0
	var scale float64 = 1.00000561024
	var falseNorthing float64 = -667.711
	var falseEasting float64 = 1500064.274

	var e2 float64 = flattening * (2.0 - flattening)
	var n float64 = flattening / (2.0 - flattening)
	var aRoof float64 = axis / (1.0 + n) * (1.0 + n*n/4.0 + n*n*n*n/64.0)
	var delta1 float64 = n/2.0 - 2.0*n*n/3.0 + 37.0*n*n*n/96.0 - n*n*n*n/360.0
	var delta2 float64 = n*n/48.0 + n*n*n/15.0 - 437.0*n*n*n*n/1440.0
	var delta3 float64 = 17.0*n*n*n/480.0 - 37*n*n*n*n/840.0
	var delta4 float64 = 4397.0 * n * n * n * n / 161280.0

	var Astar float64 = e2 + e2*e2 + e2*e2*e2 + e2*e2*e2*e2
	var Bstar float64 = -(7.0*e2*e2 + 17.0*e2*e2*e2 + 30.0*e2*e2*e2*e2) / 6.0
	var Cstar float64 = (224.0*e2*e2*e2 + 889.0*e2*e2*e2*e2) / 120.0
	var Dstar float64 = -(4279.0 * e2 * e2 * e2 * e2) / 1260.0

	var degToRad float64 = math.Pi / 180
	var lambdaZero float64 = centralMeridian * degToRad
	var xi float64 = (x - falseNorthing) / (scale * aRoof)
	var eta float64 = (y - falseEasting) / (scale * aRoof)
	var xiPrim float64 = xi -
		delta1*math.Sin(2.0*xi)*math.Cosh(2.0*eta) -
		delta2*math.Sin(4.0*xi)*math.Cosh(4.0*eta) -
		delta3*math.Sin(6.0*xi)*math.Cosh(6.0*eta) -
		delta4*math.Sin(8.0*xi)*math.Cosh(8.0*eta)
	var etaPrim float64 = eta -
		delta1*math.Cos(2.0*xi)*math.Sinh(2.0*eta) -
		delta2*math.Cos(4.0*xi)*math.Sinh(4.0*eta) -
		delta3*math.Cos(6.0*xi)*math.Sinh(6.0*eta) -
		delta4*math.Cos(8.0*xi)*math.Sinh(8.0*eta)
	var phiStar float64 = math.Asin(math.Sin(xiPrim) / math.Cosh(etaPrim))
	var deltaLambda float64 = math.Atan(math.Sinh(etaPrim) / math.Cos(xiPrim))
	var lonRadian float64 = lambdaZero + deltaLambda
	var latRadian float64 = phiStar + math.Sin(phiStar)*math.Cos(phiStar)*
		(Astar+
			Bstar*math.Pow(math.Sin(phiStar), 2)+
			Cstar*math.Pow(math.Sin(phiStar), 4)+
			Dstar*math.Pow(math.Sin(phiStar), 6))
	var lat float64 = latRadian * 180.0 / math.Pi
	var long float64 = lonRadian * 180.0 / math.Pi

	return lat, long
}
