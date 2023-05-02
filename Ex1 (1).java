package semster_b_ex1;

import java.util.Scanner;

public class ex1 {
	/**
	 * ID: 314990938
	 * 
	 * 
	 * sources:
	 * for GCD: https://www.freecodecamp.org/news/euclidian-gcd-algorithm-greatest-common-divisor/
	 * for MCD: https://www.tau.ac.il/~tsirel/dump/Static/knowino.org/wiki/Least_common_multiple.html
	 */

	public static void main(String[] args) {
		Scanner sc= new Scanner(System.in);// define scanner

		System.out.println("Enter 3 numbers for MCD:");
		// get 3 numbers
		long a=sc.nextLong();
		long b=sc.nextLong();
		long c=sc.nextLong();
		sc.close();// close scanner

		//big numbers: 9699690,30808063,87215586
		long start= System.nanoTime(); // time before running MCD algo
		long ans=(MCD(a,b,c)); // send the 3 numbers to the MCD function
		long end=System.nanoTime(); // time after running MCD algo
		System.out.println("The MCD is:"+ans);
		System.out.println("The run time is: "+(end-start)/1000 + " micro secondes");
	}
	// this function compute the MCD of 2 numbers
	public static long MCD_ab (long a, long b) {
		return a*b/GCD(a,b);
	}
	//this function compute the MCD of 3 numbers
	public static long MCD (long a, long b, long c) {
		return MCD_ab(a,MCD_ab(b,c));
	}
	// this function compute the greatest common divisor of 2 numbers
	public static long GCD (long a,long b) {
		if(b == 0)
		{
			return a;
		}
		return GCD(b, a % b);
	}
}



