/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.ij.commons;

/**
 *
 * @author airen
 */
 public class Geometry
        {
            public double shiftx, shifty, psiangle;
            
            public Geometry(double shiftx, double shifty, double psiangle)
            {
                this.shiftx = shiftx;
                this.shifty = shifty;
                this.psiangle = psiangle;
            }
        }
        