package rt90

import (
	"math"
	"testing"
)

func nearlyEqual(a, b, epsilon float64) bool {
	absA := math.Abs(a)
	absB := math.Abs(b)
	diff := math.Abs(a - b)

	if a == b {
		return true
	} else if (a == 0 || b == 0) || diff < math.SmallestNonzeroFloat64 {
		return diff < (epsilon * math.SmallestNonzeroFloat64)
	} else {
		return diff/(absA+absB) < epsilon
	}
}

func TestWGS84FromRT90(t *testing.T) {
	var tt = []struct {
		x, y, wantLat, wantLong float64
	}{
		{
			x:        6791723,
			y:        1405053,
			wantLat:  61.229502,
			wantLong: 14.037397,
		},
		{
			x:        7118097,
			y:        1575237,
			wantLat:  64.161091,
			wantLong: 17.351571,
		},
	}

	epsilon := 0.000001

	for _, test := range tt {
		lat, long := WGS84FromRT90(test.x, test.y)
		if !nearlyEqual(test.wantLat, lat, epsilon) {
			t.Errorf("unexpected latitude: got = %v, want = %v", lat, test.wantLat)
		}
		if !nearlyEqual(test.wantLong, long, epsilon) {
			t.Errorf("unexpected longitude: got = %v, want = %v", long, test.wantLong)
		}
	}
}
