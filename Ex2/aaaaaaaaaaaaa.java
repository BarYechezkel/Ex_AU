package Ex2;

import java.util.Arrays;

public class aaaaaaaaaaaaa {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		boolean ans = true;
		double EPS=0.001;
		double [] p1= {2,8,4};
		double [] p2= {2,8,4};
		int ans1 = 0;
		double count=0;
		String d= "ab ba da";
		System.out.println(d.contains("b"));
		String [] abc= d.split(" ");
		System.out.println( Arrays.toString(abc));
		for (int i=0; i<abc.length; i++) {
			System.out.println(abc[i]);
		}
		if (p1.length == p2.length)
		{
			for (int i=0; i< p1.length;i++)
			{
				double a=Math.abs(p1[i]-p2[i]);
				if (a<=EPS) {
					ans1++;
				}
			}
			if (ans1 == p1.length) {
				ans= true;
			}
		}
System.out.println(ans);
	}

}
