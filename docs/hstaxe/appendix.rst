.. _appendix:

========
Appendix
========


Virtual slits
=============

The slit length *sl*, the slit width *sw* and the
orientation *so* of the *virtual slit*, optimized to maintain the
spectral resolution in slitless spectroscopy, is computed from the major
axis size a_img, the minor axis size b_img and the major axis angle
theta_img for each object (the three quantities are given in the
SExtractor catalogue columns A_IMAGE, B_IMAGE and THETA_IMAGE).

Without loss of generality we assume that the dispersion direction is
parallel to the x-axis, which means the angle between the dispersion
direction and the major axis angle (defined with respect to the x-axis)
is theta_img. It is then:

.. math::

   \begin{aligned}
   A_{11}     & = & (\cos({theta\_img}) / a\_img)^2 + (\sin{(theta\_img)} / b\_img)^2\\
   A12        & = & \cos{(theta\_img)} * \sin{(theta\_img)} * (1.0/a\_img^2 - 1.0/b\_img^2)\\
   A_{22}     & = & (\sin({theta\_img}) / a\_img)^2 + (\cos{(theta\_img)} / b\_img)^2\\
   \alpha     & = & \arctan{(A12 / A11)}\\\nonumber \\
   sl         & = &  \sqrt{A11} * a\_img * b\_img / \cos{(\alpha)} \\
   sw         & = &  1.0 / \sqrt{A11} \\
   so         & = &  \alpha + 90.0\end{aligned}

For further details see [FREUDLING]_

