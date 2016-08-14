#include "../Geometry_util/Geometry_util.h"

#include <vector>

using namespace std;


void fill_up_fragment_torsion_angles ( 
 const double * coord_array,
 const int length, 
 vector < double > & torsion_set, 
 const char radian_or_gradus_flag  )
{
	torsion_set.clear();

	double p1[3],p2[3],p3[3],p4[3];

	for (int ii=2;ii<3*length-1; ii++)
	{
	   	
		p1 [0] = coord_array [3*(ii-2)		];
		p1 [1] = coord_array [3*(ii-2) + 1	];
		p1 [2] = coord_array [3*(ii-2) + 2	];


		p2 [0] = coord_array [3*(ii-1)		];
		p2 [1] = coord_array [3*(ii-1) + 1	];
		p2 [2] = coord_array [3*(ii-1) + 2	];


		p3[0] = coord_array [3*ii			];
		p3[1] = coord_array [3*ii		+ 1	];
		p3[2] = coord_array [3*ii		+ 2	];
		
		p4 [0] = coord_array [3*(ii+1)		];
		p4 [1] = coord_array [3*(ii+1) + 1	];
		p4 [2] = coord_array [3*(ii+1) + 2	];


		torsion_set .push_back(  Geometry_util::diheral_by_four_poins ( p1,p2,p3,p4 )  );
	}

	if ( radian_or_gradus_flag == 'd' || radian_or_gradus_flag == 'd' )
		for (int ii=0;ii <torsion_set.size() ; ii++)
			torsion_set[ii] *= 180.0/ Geometry_util::Pythagorean_Number ();

}
