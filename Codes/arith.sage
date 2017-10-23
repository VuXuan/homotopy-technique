def exp(self, prec=None):
        r"""
        Return exp of this power series to the indicated precision.

        INPUT:


        -  ``prec`` - integer; default is
           ``self.parent().default_prec``


        ALGORITHM: See :meth:`.solve_linear_de`.

        .. note::

           - Screwy things can happen if the coefficient ring is not a
             field of characteristic zero. See :meth:`.solve_linear_de`.

        AUTHORS:

        - David Harvey (2006-09-08): rewrote to use simplest possible "lazy" algorithm.

        - David Harvey (2006-09-10): rewrote to use divide-and-conquer
          strategy.

        - David Harvey (2006-09-11): factored functionality out to
          solve_linear_de().

        - Sourav Sen Gupta, David Harvey (2008-11): handle constant
          term

        EXAMPLES::

            sage: R.&lt;t&gt; = PowerSeriesRing(QQ, default_prec=10)

        Check that `\exp(t)` is, well, `\exp(t)`::

            sage: (t + O(t^10)).exp()
            1 + t + 1/2*t^2 + 1/6*t^3 + 1/24*t^4 + 1/120*t^5 + 1/720*t^6 + 1/5040*t^7 + 1/40320*t^8 + 1/362880*t^9 + O(t^10)

        Check that `\exp(\log(1+t))` is `1+t`::

            sage: (sum([-(-t)^n/n for n in range(1, 10)]) + O(t^10)).exp()
            1 + t + O(t^10)

        Check that `\exp(2t + t^2 - t^5)` is whatever it is::

            sage: (2*t + t^2 - t^5 + O(t^10)).exp()
            1 + 2*t + 3*t^2 + 10/3*t^3 + 19/6*t^4 + 8/5*t^5 - 7/90*t^6 - 538/315*t^7 - 425/168*t^8 - 30629/11340*t^9 + O(t^10)

        Check requesting lower precision::

            sage: (t + t^2 - t^5 + O(t^10)).exp(5)
            1 + t + 3/2*t^2 + 7/6*t^3 + 25/24*t^4 + O(t^5)

        Can't get more precision than the input::

            sage: (t + t^2 + O(t^3)).exp(10)
            1 + t + 3/2*t^2 + O(t^3)

        Check some boundary cases::

            sage: (t + O(t^2)).exp(1)
            1 + O(t)
            sage: (t + O(t^2)).exp(0)
            O(t^0)

        Handle nonzero constant term (fixes :trac:`4477`)::

            sage: R.&lt;x&gt; = PowerSeriesRing(RR)
            sage: (1 + x + x^2 + O(x^3)).exp()
            2.71828182845905 + 2.71828182845905*x + 4.07742274268857*x^2 + O(x^3)

        ::

            sage: R.&lt;x&gt; = PowerSeriesRing(ZZ)
            sage: (1 + x + O(x^2)).exp()
            Traceback (most recent call last):
            ...
            ArithmeticError: exponential of constant term does not belong to coefficient ring (consider working in a larger ring)

        ::

            sage: R.&lt;x&gt; = PowerSeriesRing(GF(5))
            sage: (1 + x + O(x^2)).exp()
            Traceback (most recent call last):
            ...
            ArithmeticError: constant term of power series does not support exponentiation
        """
        if prec is None:
            prec = self._parent.default_prec()

        t = self.derivative().solve_linear_de(prec)

        if not self[0].is_zero():
            try:
                C = self[0].exp()
            except AttributeError:
                raise ArithmeticError("constant term of power series does not support exponentiation")

            if C.parent() is not self.base_ring():
                raise ArithmeticError("exponential of constant term does not belong to coefficient ring (consider working in a larger ring)")

            t = C * t

        return t
