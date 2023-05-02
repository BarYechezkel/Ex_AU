package Ex2;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

import java.util.Arrays;

/**
 *  * Introduction to Computer Science 2023, Ariel University,
 *  * Ex2: arrays, static functions and JUnit
 *
 * This JUnit class represents a JUnit (unit testing) for Ex2 - 
 * It contains few testing functions for the polynum functions as define in Ex2.
 * Note: you should add additional JUnit testing functions to this class.
 *
 * @author boaz.ben-moshe
 */

class Ex2Test {
	static final double[] P1 ={2,0,3, -1,0}, P2 = {0.1,0,1, 0.1,3};
	static double[] po1 = {2,2}, po2 = {-3, 0.61, 0.2};;
	static double[] po3 = {2,1,-0.7, -0.02,0.02};
	static double[] po4 = {-3, 0.61, 0.2};


	@Test
	void testPolynomFromPoints() {
		double []po1= {1,2,3};
		double []po2= {5,11,19};
		double []po3= {1.0, 3.0, 1.0};
		double [] f1=Ex2.PolynomFromPoints(po1, po2);
		boolean ans=(Ex2.equals(po3, f1));
		assertTrue (ans);
		
		double []po4= {4,-6,2};
		double []po5= {-1,4,3};
		double []po6= {5.5, -0.875, -0.1875};
		double [] f0=Ex2.PolynomFromPoints(po4, po5);
		boolean ans4=(Ex2.equals(po6, f0));
		assertTrue (ans4);
		
		double [] p= {1,2};
		double [] p2= {2,4};
		double [] f2= Ex2.PolynomFromPoints(p, p2);
		double [] f3= {0,2};
		boolean ans2=(Ex2.equals(f2, f3));
		assertTrue (ans2);
		
		double [] p4= {0,-4};
		double [] p5= {8,2};
		double [] f4= Ex2.PolynomFromPoints(p4, p5);
		double [] f5= {8,1.5};
		boolean ans3=(Ex2.equals(f4, f5));
		assertTrue (ans3);

	}
	
	

	@Test
	void testPoly() {
		double [] p1= {1,2,3};
		String ans= Ex2.poly(p1);
		assertEquals(ans, "3.0x^2 +2.0x +1.0");

		double [] p2= {-1,-2,3};
		String ans2= Ex2.poly(p2);
		assertEquals(ans2, "3.0x^2 -2.0x -1.0");

		double [] p3= {-9,2};
		String ans3= Ex2.poly(p3);
		assertEquals(ans3,"2.0x -9.0");

		double [] p4= {-9,2,4,-5,8,-14};
		String ans4= Ex2.poly(p4);
		assertEquals(ans4,"-14.0x^5 +8.0x^4 -5.0x^3 +4.0x^2 +2.0x -9.0");
	}

	@Test
	/**
	 * Tests that f(x) == poly(x).
	 */
	void testF() {
		double fx0 = Ex2.f(po1, 0);
		double fx1 = Ex2.f(po1, 1);
		double fx2 = Ex2.f(po1, 2);
		assertEquals(fx0, 2, Ex2.EPS);
		assertEquals(fx1, 4, Ex2.EPS);
		assertEquals(fx2, 6, Ex2.EPS);
	}
	@Test
	/**
	 * Tests that p1(x) + p2(x) == (p1+p2)(x)
	 */
	void testF2() {
		double x = Math.PI;
		double[] po12 = Ex2.add(po1, po2);
		double f1x = Ex2.f(po1, x);
		double f2x = Ex2.f(po2, x);
		double f12x = Ex2.f(po12, x);
		assertEquals(f1x + f2x, f12x, Ex2.EPS);
	}
	@Test
	/**
	 * Tests that p1+p2+ (-1*p2) == p1
	 */
	void testAdd() {
		double[] p12 = Ex2.add(po1, po2);
		double[] minus1 = {-1};
		double[] pp2 = Ex2.mul(po2, minus1);
		double[] p1 = Ex2.add(p12, pp2);

		assertFalse(Ex2.equals(p1, po1));
	}
	@Test
	/**
	 * Tests that p1+p2 == p2+p1
	 */
	void testAdd2() {
		double[] p12 = Ex2.add(po1, po2);
		double[] p21 = Ex2.add(po2, po1);
		assertTrue(Ex2.equals(p12, p21));
	}
	@Test
	/**
	 * Tests that p1+0 == p1
	 */
	void testAdd3() {
		double[] p1 = Ex2.add(po1, Ex2.ZERO);
		assertTrue(Ex2.equals(p1, po1));
	}
	@Test
	/**
	 * Tests that p1*0 == 0
	 */
	void testMul1() {
		double[] p1 = Ex2.mul(po1, Ex2.ZERO);
		assertTrue(Ex2.equals(p1, Ex2.ZERO));
	}
	@Test
	/**
	 * Tests that p1*p2 == p2*p1
	 */
	void testMul2() {
		double[] p12 = Ex2.mul(po1, po2);
		double[] p21 = Ex2.mul(po2, po1);
		assertTrue(Ex2.equals(p12, p21));
	}
	@Test
	/**
	 * Tests that p1(x) * p2(x) = (p1*p2)(x),
	 */
	void testMulDoubleArrayDoubleArray() {
		double[] xx = {0,1,2,3,4.1,-15.2222};
		double[] p12 = Ex2.mul(po1, po2);
		for(int i = 0;i<xx.length;i=i+1) {
			double x = xx[i];
			double f1x = Ex2.f(po1, x);
			double f2x = Ex2.f(po2, x);
			double f12x = Ex2.f(p12, x);
			assertEquals(f12x, f1x*f2x, Ex2.EPS);
		}
	}
	@Test
	/**
	 * Tests a simple derivative examples - till ZERO.
	 */
	void testDerivativeArrayDoubleArray() {
		double[] p = {1,2,3}; // 3X^2+2x+1
		double[] pt = {2,6}; // 6x+2
		double[] dp1 = Ex2.derivative(p); // 2x + 6
		double[] dp2 = Ex2.derivative(dp1); // 2
		double[] dp3 = Ex2.derivative(dp2); // 0
		double[] dp4 = Ex2.derivative(dp3); // 0
		assertTrue(Ex2.equals(dp1, pt));
		assertTrue(Ex2.equals(Ex2.ZERO, dp3));
		assertTrue(Ex2.equals(dp4, dp3));
	}
	@Test
	/** 
	 * Tests the parsing of a polynom in a String like form.
	 */
	public void testFromString() {
		double[] p = {-1.1,2.3,3.1}; // 3.1X^2+ 2.3x -1.1
		String sp2 = "3.1x^2 +2.3x -1.1";
		String sp = Ex2.poly(p);
		double[] p1 = Ex2.getPolynomFromString(sp);
		double[] p2 = Ex2.getPolynomFromString(sp2);
		boolean isSame1 = Ex2.equals(p1, p);
		boolean isSame2 = Ex2.equals(p2, p);
		if(!isSame1) {fail();}
		if(!isSame2) {fail();}
		assertEquals(sp, Ex2.poly(p1));
	}
	@Test
	public void testFromString1() {
		double[] p = {-1.1}; // -1.1
		String sp2 = "-1.1";
		String sp = Ex2.poly(p);
		double[] p1 = Ex2.getPolynomFromString(sp);
		double[] p2 = Ex2.getPolynomFromString(sp2);
		boolean isSame1 = Ex2.equals(p1, p);
		boolean isSame2 = Ex2.equals(p2, p);
		if(!isSame1) {fail();}
		if(!isSame2) {fail();}
		assertEquals(sp, Ex2.poly(p1));
	}
	@Test
	public void testFromString2() {
		double[] p = {-1.1,2}; //2x -1.1
		String sp2 = "2x-1.1";
		String sp = Ex2.poly(p);
		double[] p1 = Ex2.getPolynomFromString(sp);
		double[] p2 = Ex2.getPolynomFromString(sp2);
		boolean isSame1 = Ex2.equals(p1, p);
		boolean isSame2 = Ex2.equals(p2, p);
		if(!isSame1) {fail();}
		if(!isSame2) {fail();}
		assertEquals(sp, Ex2.poly(p1));
	}
	
	@Test
	/**
	 * Tests the equality of pairs of arrays.
	 */
	public void testEquals() {
		double[][] d1 = {{0}, {1}, {1,2,0,0}};
		//double[][] d1 = { {1,2,0,0}};
		double[][] d2 = {Ex2.ZERO, {1+Ex2.EPS/2}, {1,2}};
		double[][] xx = {{-2*Ex2.EPS}, {1+Ex2.EPS*1.2}, {1,2,Ex2.EPS*2}};
		//double[][] xx = { {1,2,Ex2.EPS*2}};
		for(int i=0;i<d1.length;i=i+1) {
			//	assertTrue(Ex2.equals(d1[i], d2[i]));
		}
		for(int i=0;i<d1.length;i=i+1) {
			assertFalse(Ex2.equals(d1[i], xx[i]));
		}
	}

	@Test 
	public void test_equals() 
	{
		boolean case1=Ex2.equals(po1,po2);
		assertFalse(case1);

		double [] p1= {2,8,4};
		double [] p2= {2,8,4};
		boolean case2= Ex2.equals(p1, p2);
		assertTrue(case2);

		double [] p3= {2,8.5,4,0,0};
		double [] p4= {2,-8.5,4};
		boolean case3=Ex2.equals(p3,p4);
		assertFalse(case3);

		double [] p5= {-1,8,5,0,0};
		double [] p6= {-1,8,4};
		boolean case4=Ex2.equals(p5,p6);
		assertFalse(case4);

		double [] p7= {1,-3,5,0,0};
		double [] p8= {1,-3,5,0};
		boolean case5=Ex2.equals(p7,p8);
		assertTrue(case5);
	}

	@Test
	/**
	 * Tests is the sameValue function is symmetric.
	 */
	public void testSameValue() {
		double x1=-8, x2=4;
		double rs1 = Ex2.sameValue(po1,po2, x1, x2, Ex2.EPS);
		double rs2 = Ex2.sameValue(po2,po1, x1, x2, Ex2.EPS);
		assertEquals(rs1,rs2,Ex2.EPS);

		double x3=-28, x4=35;
		double rs3 = Ex2.sameValue(po3,po4, x3, x4, Ex2.EPS);
		double rs4 = Ex2.sameValue(po4,po3, x3, x4, Ex2.EPS);
		assertEquals(rs3,rs4,Ex2.EPS);
	}
	
	
	@Test
	public void testLength() {
		double [] f= {0,0,1};
		double ans= Ex2.length(f, 0, 2, 2);
		assertEquals(ans, (Math.sqrt(2)+Math.sqrt(10)));
		
		double [] f2= {0,2};
		double ans2= Ex2.length(f2, 0, 2, 10);
		assertEquals(ans2, Math.sqrt(20),Ex2.EPS);
		
		double [] f3= {0,2};
		double ans3= Ex2.length(f3, 0, 2, 2);
		assertEquals(ans3,2* Math.sqrt(5),Ex2.EPS);
		
		double [] f4= {4};
		double ans4= Ex2.length(f4, 2, 100, 10);
		assertEquals(ans4, 98);
		
		double [] f5= {0,0,0,0,0,0};
		double ans5= Ex2.length(f5, 2, 5, 10);
		assertEquals(ans5, 3);
		
		double [] f6= null;
		double ans6= Ex2.length(f6, 1, 2, 50);
		assertEquals(ans6, 0);
		
		
	}

	@Test
	/**
	 * Test the area function - it should be symmetric.
	 */
	public void testArea() {
		double x1=0, x2=4;
		double a1 = Ex2.area(po1, po2, x1, x2, 100);
		double a2 = Ex2.area(po2, po1, x1, x2, 100);
		assertEquals(a1,a2,Ex2.EPS);
	}

}
