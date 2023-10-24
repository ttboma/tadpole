#ifndef _SIMULATED_ANEALING_H
#define _SIMULATED_ANEALING_H

namespace syc {

	namespace SA::v_01 {

		template <
			typename Sol, 
			typename T,
			typename SchulePolicy1,
			typename SchulePolicy2,
			typename SchulePolicy3
		>
		Sol 
		simulated_annealing( 
			Sol& sNow, T currT, T const endingT, 
			SchulePolicy1	thermal_equilibrium, 
			SchulePolicy2	decrease,
			SchulePolicy3	accept
		) {
			Sol sBest = sNow;
			while ( currT > endingT ) 
			{
				while ( !thermal_equilibrium( currT ) ) 
				{
					auto sNext  = sNow;
				 	sNext.perturb();

					auto deltaC = sNext.cost() - sNow.cost(); 	
					if (deltaC < 0) 
					{
						sNow = sNext;
						if (sNext.cost() < sBest.cost()) 
						{
							sBest = sNext;
						}
					}
					else if (accept( currT, deltaC )) 
					{
						sNow = sNext;
					}
				}
				decrease( currT );
			}
			return sBest;
		}
	}			
}

#endif
