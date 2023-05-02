package Ex2;

import java.util.Arrays;

/**
 * Introduction to Computer Science 2023, Ariel University,
 * Ex2: arrays, static functions and JUnit
 *
 * This class represents a set of functions on a polynom - represented as array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynom: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code here")
 *
 * @author boaz.benmoshe
 */
public class Ex2 {
	/** Epsilon value for numerical computation, it serves as a "close enough" threshold. */
	public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
	/** The zero polynom is represented as an array with a single (0) entry. */
	public static final double[] ZERO = {0};
	/**
	 * Computes the f(x) value of the polynom at x.
	 * @param poly
	 * @param x
	 * @return f(x) - the polynom value at x.
	 */
	public static double f(double[] poly, double x) {
		double ans = 0;
		for(int i=0;i<poly.length;i++) {
			double c = Math.pow(x, i);
			ans +=c*poly[i];
		}
		return ans;
	}
	/** Given a polynom (p), a range [x1,x2] and an epsilon eps. 
	 * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps, 
	 * assuming p(x1)*p(x1) <= 0. 
	 * This function should be implemented recursively.
	 * @param p - the polynom
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
	 */
	public static double root_rec(double[] p, double x1, double x2, double eps) {
		double f1 = f(p,x1);
		double f2 = f(p,x2);
		double x12 = (x1+x2)/2;
		double f12 = f(p,x12);
		if (f1*f2<=0 && Math.abs(f12)<eps) {return x12;}
		if(f12*f1<=0) {return root_rec(p, x1, x12, eps);}
		else {return root_rec(p, x12, x2, eps);}
	}
	/**
	 * This function computes a polynomial representation from a set of 2D points on the polynom.
	 * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
	 * Note: this function only works for a set of points containing up to 3 points, else returns null.
	 * @param xx
	 * @param yy
	 * @return an array of doubles representing the coefficients of the polynom.
	 */
	public static double[] PolynomFromPoints(double[] xx, double[] yy) {
		double [] ans = null;
		int lx = xx.length;
		int ly = yy.length;
		if(xx!=null && yy!=null && lx==ly && lx>1 && lx<4) {
			if (xx[0]==xx[1]) {//m= (y2-y1)/(x2-x1); x2-x1 cannot be '0'.
				return null; //
			}
			// if the polynom have 2 monoms
			if (lx==ly && lx==2) {
				ans= new double [2];
				// y=mx+n
				double x1= xx[0], x2= xx[1];
				double y1= yy[0], y2= yy[1];
				double m=0,n=0;
				//compute the m and n
				m= (y2-y1)/(x2-x1);
				n=y1-(m*x1);
				//put the a,b,c in the appropriate place in the array
				ans [0]= n;
				ans [1]= m;
			}
			else{// if the polynom have 3 monoms
				ans= new double [3];
				double x1= xx[0],x2= xx[1],x3= xx[2];
				double y1= yy[0],y2= yy[1],y3= yy[2];
				double a=0,b=0,c=0;
				if (((x1-x2)*(x1-x3)*(x2-x3))==0) {
					return null; // divide by zero
				}
				//y=ax^2+bx+c
				a= ((x1*(y3-y2)+x2*(y1-y3)+x3*(y2-y1))/((x1-x2)*(x1-x3)*(x2-x3)));
				b= (((y2-y1)/(x2-x1))-(a*(x1+x2)));
				c= ((y1-(a*(x1*x1))-(b*x1)));
				//put the a,b,c in the appropriate place in the array
				ans[0]=c;
				ans[1]=b;
				ans[2]=a;
			}

			////////////////////
		}
		return ans;
	}
	/** Two polynoms are equal if and only if the have the same values f(x) for 1+n values of x, 
	 * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
	 * @param p1 first polynom
	 * @param p2 second polynom
	 * @return true iff p1 represents the same polynom as p2.
	 */
	public static boolean equals(double[] p1, double[] p2) {
		boolean ans = true;
		int ans1 = 0;
		double count=0;
		ans= false;
		// if one of the polynoms or both of them are null
		if (p1==null || p2==null) {
			ans = false;
		}
		// check if the polynoms have the same size 
		if (p1.length == p2.length)
		{
			for (int i=0; i< p1.length;i++)
			{
				// compare the the arrays elements of the two polynoms up to EPS
				if ((double)Math.abs(p1[i]-p2[i])<EPS) {
					ans1++;
				}
			}
			if (ans1 == p1.length) {
				ans= true;
			}
		}
		// if p1 bigger than p2
		if (p1.length < p2.length) 
		{
			for (int i=0; i< p1.length;i++)
			{
				//check if the subtraction between every element in the array is less then EPS, if its true ans1 increased by one
				if ((Math.abs(p1[i]-p2[i]))<EPS) {
					ans1++;
				}
			}
			int tmp = p2.length- p1.length;
			count=0;
			// count the elements in p2, from the p1 size to the end of p2
			for (int j=0; j < tmp ; j++)
			{
				count= count + p2[p1.length+j];

			}
			// if count smaller than EPS and the counter ans1 equals to p1 size- the polynoms are equals.
			if ((count < EPS) && (ans1== p1.length)) {
				ans= true;
			}

		}
		if (p1.length > p2.length) // if p2 bigger than p1
		{
			for (int i=0; i< p2.length;i++)
			{
				//check if the subtraction between every element in the array is less then EPS, if its true ans1 increased by one
				if ((Math.abs(p1[i]-p2[i]))<EPS)
				{
					ans1++;
				}
			}

			int tmp = p1.length - p2.length;
			count=0;
			// count the elements in p1, from the p2 size to the end of p1
			for (int j=0; j< tmp ; j++)
			{
				count= count + p1[p2.length+j];

			}
			//if count smaller than EPS and the counter ans1 equals to p2 size- the polynoms are equals.
			if ((count < EPS) && (ans1== p2.length)) { 
				ans= true;
			}
		}

		////////////////////
		return ans;
	}

	/** 
	 * Computes a String representing the polynom.
	 * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
	 * @param poly the polynom represented as an array of doubles
	 * @return String representing the polynom: 
	 */
	public static String poly(double[] poly) {
		String ans = "";
		if(poly.length==0) {ans="0";}
		else {
			// add you code here
			String y = "";// 'hezka' for monom> 2;
			String[] temp= new String[poly.length];// array of strings for the given polynom
			if (poly.length==1)// if poly is one number
			{
				temp[0]=Double.toString(poly[0]);
			}
			if (poly.length>1)
			{
				if (poly[0]!=0)
				{
					if (poly[0]<0)
					{
						temp[0]=Double.toString(poly[0]);
					}
					if (poly[0]>0)
					{
						temp[0]="+" +Double.toString(poly[0]);
					}

				}
				if (poly[1]!=0)//if poly[1] is not '0'
					temp[1]=Double.toString(poly[1])+"x";
				if (poly.length>2)// if the numbers near poly[1] (poly[2]) is '0'.
				{
					if (poly[1]<0)
						temp[1]=Double.toString(poly[1])+"x";
					if (poly[1]>0)
						temp[1]="+" +Double.toString(poly[1])+"x";
				}
				if ((poly[poly.length-1]!=0)&&(poly.length>2))// handle the case that the number with the bigeest 'hezka' is positive
				{
					y=Integer.toString(poly.length-1);
					temp[poly.length-1]=Double.toString(poly[poly.length-1])+"x^"+y;

				}

				for (int i=2;i<poly.length-1;i++)// put the monoms of hezka>2 in the same location in temp
				{
					if(poly[i]>0)
					{
						y=Integer.toString(i);
						temp[i]="+"+Double.toString(poly[i])+"x^"+y;
					}
					if(poly[i]<0)
					{
						y=Integer.toString(i);
						temp[i]=Double.toString(poly[i])+"x^"+y;	
					}	
				}
			}
			for (int i=poly.length-1; i>0; i--)// sum all the monoms that not '0' in temp array.
			{
				if(poly[i]!=0)
					ans+= temp[i]+" ";
			}
			ans=ans+temp[0];

			////////////////////
		}
		return ans;
	}
	/**
	 * Given two polynoms (p1,p2), a range [x1,x2] and an epsilon eps. This function computes an x value (x1<=x<=x2)
	 * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
	 * @param p1 - first polynom
	 * @param p2 - second polynom
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
	 */
	public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
		double ans = x1;
		// find the x value of the intersection point between 2 functions
		double x12 = (x1+x2)/2;
		while (Math.abs((f(p1, x12))-(f(p2, x12)))> eps)
		{
			if ((f(p1, x12) - f(p2, x12))*((f(p1, x1) -f(p2, x1)))>=0)
			{
				x1=x12;
			}
			else 
			{
				x2=x12;
			}	
			x12=(x1+x2)/2;
		}
		ans= x12;
		return ans;
	}
	/**
	 * Given a polynom (p), a range [x1,x2] and an integer with the number (n) of sample points. 
	 * This function computes an approximation of the length of the function between f(x1) and f(x2) 
	 * using n inner sample points and computing the segment-path between them.
	 * assuming x1 < x2. 
	 * This function should be implemented iteratively (none recursive).
	 * @param p - the polynom
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfSegments - (A positive integer value (1,2,...).
	 * @return the length approximation of the function between f(x1) and f(x2).
	 */
	public static double length(double[] p, double x1, double x2, int numberOfSegments) {
		double ans = x1;
		// add you code here
		double sum=0;
		//check if the polynom is empty
		if (p==null) {
			return 0;
		}
		//compute the size of every segment
		double segmentSize=Math.abs(x2-x1)/numberOfSegments;

		for(int i=0; i < numberOfSegments ;i++) {
			//distance: sqrt((y2-y1)^2+(x2-x1)^2)
			double xx2=x1+((i+1)*segmentSize);
			double xx1=x1+(i*segmentSize);
			double y2=f(p,xx2);
			double y1=f(p,xx1);

			double a= Math.pow(y2-y1,2);
			double b= Math.pow(xx2-xx1,2);
			//sum all the size of the segments 
			sum=sum+Math.sqrt(a+b);
		}
		////////////////////
		return sum;
	}

	/**
	 * Given two polynoms (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids between the functions (number of samples in on each polynom). 
	 * This function computes an approximation of the area between the polynoms within the x-range.
	 * The area is computed using Riemann's like integral (https://en.wikipedia.org/wiki/Riemann_integral)
	 * @param p1 - first polynom
	 * @param p2 - second polynom
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
	 * @return the approximated area between the two polynoms within the [x1,x2] range.
	 */
	public static double area(double[] p1,double[]p2, double x1, double x2, int numberOfTrapezoid) {
		double ans = 0;
		double sum =0;
		// create arrays of the x's that represent the x's value of each trapz heights
		double [] xx= new double [numberOfTrapezoid+1];
		// the size of the height of each trapez
		double nSize= (x2-x1)/numberOfTrapezoid;
		//put the x's value in the xx array
		for (int i=0; i<xx.length; i++) {
			xx[i]=x1+nSize*i;
		}
		for (int i=1; i<=numberOfTrapezoid ;i++) {
			//check if there is intersection point between 2 x's, if it is, compute area of 2 triangles
			if((f(p2,xx[i])-f(p1,xx[i]))*(f(p2,xx[i-1])-f(p1,xx[i-1]))<0){
				//compute the intersection point
				double x12=sameValue(p1, p2, xx[i-1], xx[i], EPS);
				double basis1=0;
				double basis2=0;
				//compute the bases of the two triangle
				if (f(p1,xx[i-1])>f(p2,xx[i-1])) {
					basis1=f(p1,xx[i-1])-f(p2,xx[i-1]);
				}
				else {
					basis1=f(p2,xx[i-1])-f(p1,xx[i-1]);
				}
				if (f(p1,xx[i])>f(p2,xx[i])) {
					basis2=f(p1,xx[i])-f(p2,xx[i]);
				}
				else {
					basis2=f(p2,xx[i])-f(p1,xx[i]);
				}
				//compute the area of the triangles
				double triangle1= ((x12-xx[i-1])*basis1)/2;
				double triangle2= ((xx[i]-x12)*basis2)/2;
				sum=sum +(triangle1+triangle2);


			}
			//check if there is no intersection point between 2 x's, compute area of trapez
			else {
				// f1(x)>f2(x)
				if((f(p2,xx[i])-f(p1,xx[i])<=EPS)){
					double heightA=f(p1,xx[i-1])-f(p2,xx[i-1]);
					double heightB=f(p1,xx[i])-f(p2,xx[i]);
					sum=sum + ((heightA+heightB)*nSize)/2;
				}
				// f1(x)<f2(x)
				else {
					double heightA=f(p2,xx[i-1])-f(p1,xx[i-1]);
					double heightB=f(p2,xx[i])-f(p1,xx[i]);
					sum=sum+ ((heightA+heightB)*nSize)/2;
				}
			}
		}
		ans=sum;
		return ans;
	}
	/**
	 * This function computes the array representation of a polynom from a String
	 * representation. Note:given a polynom represented as a double array,  
	 * getPolynomFromString(poly(p)) should return an array equals to p.
	 * 
	 * @param p - a String representing polynom.
	 * @return
	 */
	public static double[] getPolynomFromString(String p) {
		double [] ans = ZERO;//  -1.0x^2 +3.0x +2.0
		// add you code here
		String p_new= p.replace(" ", "");
		p_new= p_new.replace("-", "+-");
		int max=1;
		int power=0;
		String[] tmp=p.split("((?>=-)|(?=-)|(?>=\\+)|(?=\\+))");//split by + & - and create array of monoms.
		int numOfMonoms=tmp.length;
		if(numOfMonoms==1) {
			max=0;
		}
		if(numOfMonoms>2) {
			p_new=p_new.substring(p_new.indexOf('^')+1);//find the 'hezka'
			max=Character.getNumericValue(p_new.charAt(0));
		}
		double[] arrc= new double[max+1];
		for (int i=0; i<numOfMonoms;i++) {
			if (tmp[i].indexOf('^')==-1 && tmp[i].indexOf("x")==-1) {
				arrc[0]=Double.parseDouble(tmp[i]);
			}
			else
				if (tmp[i].indexOf('^')==-1) {
					String[] splitByX= tmp[i].split("x");
					arrc[1]=Double.parseDouble(splitByX[0]);
				}
				else {
					String[] splitByX= tmp[i].split("x");
					splitByX[1]=splitByX[1].substring((splitByX[1].indexOf('^')+1));
					String num= splitByX[1];
					power=(int)Double.parseDouble(num);
					arrc[power]=Double.parseDouble(splitByX[0]);;
				}
		}
		return arrc;
	}
	/**
	 * This function computes the polynom which is the sum of two polynoms (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static double[] add(double[] p1, double[] p2) {
		double [] ans = ZERO;//
		// add you code here	
		if (p1.length>p2.length) {
			ans= new double [p1.length];
			for (int i=0; i<p2.length; i++) {
				ans[i]=p1[i]+p2[i];
			}
			for (int i=p2.length; i<p1.length; i++) {
				ans[i]= p1[i];
			}
		}
		if (p1.length<p2.length) {
			ans= new double [p2.length];
			for (int i=0; i<p1.length; i++) {
				ans[i]=p1[i]+p2[i];
			}
			for (int i=p1.length; i<p2.length; i++) {
				ans[i]= p2[i];
			}
		}
		////////////////////
		return ans;
	}
	/**
	 * This function computes the polynom which is the multiplication of two polynoms (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static double[] mul(double[] p1, double[] p2) {
		double [] ans = ZERO;//
		// add you code here
		ans= new double [(p2.length-1)+(p1.length-1)+1];
		for (int i=0; i<p1.length; i++) {
			for (int j=0; j<p2.length; j++) {
				ans[i+j]+=p1[i]*p2[j];
			}
		}
		////////////////////
		return ans;
	}
	/**
	 * This function computes the derivative polynom:.
	 * @param po
	 * @return
	 */
	public static double[] derivative (double[] po) {
		double [] ans = ZERO;//
		if (po.length<=1) {
			return ans;
		}
		else {
			ans= new double [po.length-1];
			for (int i=0; i<ans.length; i++) {
				ans[i]=po[i+1]*(i+1);
			}
		}
		return ans;
	}
}
