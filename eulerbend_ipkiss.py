""" Generate an euler bend with IPKISS """

from __future__ import division

import numpy as np
from scipy.special import fresnel

import technologies.silicon_photonics
import ipkiss3.all as i3


def _euler_bend_center_line_shape(R=10, theta=np.pi / 2, num_points=1000):
    """
    Args:
        R (float): minimal bend radius
        theta (float): final angle
        num_points (int): resolution of the shape
    """
    L = 2 * R * theta / 2  # HALF of total length
    s = np.linspace(0, L, num_points // 2)
    f = np.sqrt(np.pi * R * L) + 1e-18  # for numerical stability
    y1, x1 = fresnel(s / f)
    # first, rotate by the final angle
    x2, y2 = np.dot(
        np.array([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]]),
        np.stack([x1, y1], 0),
    )
    # then, flip along the x-axis (and reverse direction of the curve):
    x2, y2 = -x2[::-1], y2[::-1]
    # then translate from (x2[0], y2[0]) to (x1[-1], y1[-1])
    x2, y2 = x2 - x2[0] + x1[-1], y2 - y2[0] + y1[-1]
    x = f * np.concatenate([x1, x2], 0)
    y = f * np.concatenate([y1, y2], 0)
    return zip(x, y)


class EulerBend(i3.PCell):
    """ An Euler Bend for IPKISS """

    _name_prefix = "EULER_BEND"

    trace_template = i3.WaveguideTemplateProperty(
        doc="trace template for the waveguide"
    )
    def _default_trace_template(self):
        return i3.TECH.PCELLS.WG.DEFAULT

    waveguide = i3.ChildCellProperty(doc="waveguide")
    def _default_waveguide(self):
        return i3.RoundedWaveguide(trace_template=self.trace_template)

    class Layout(i3.LayoutView):
        """ Layout for the IPKISS Euler Bend """

        min_radius = i3.PositiveNumberProperty(doc="minimum Euler bend radius")
        def _default_min_radius(self):
            if self.end_point is None:
                return 10
            x, y = self.end_point
            theta = 2 * np.arctan2(
                y, x
            )  # the final angle is twice the angle between start and finish.
            x1, y1 = euler_bend(R=1, theta=theta, num_points=4)
            return x / x1[-1]

        end_angle = i3.PositiveNumberProperty(doc="angle at the end of the Euler bend.")
        def _default_end_angle(self):
            if self.end_point is None:
                return 90
            x, y = self.end_point
            return 2 * np.arctan2(y, x) * i3.RAD2DEG

        num_points = i3.PositiveNumberProperty(
            doc="resolution of the euler bend shape.",
            default=1000,
        )

        end_point = i3.Coord2Property(
            doc=(
                "(optional) coordinates of the end of the Euler bend.\n"
                "if given, this will override `min_radius` and `end_angle`."
            ),
            default=None,
        )

        def _default_trace_template(self):
            trace = self.cell.trace_template.get_default_view(i3.LayoutView)
            trace.set(
                core_width=1,
                cladding_width=5,
            )
            return trace

        def _default_waveguide(self):
            waveguide = self.cell.waveguide.get_default_view(i3.LayoutView)
            center_line = _euler_bend_center_line_shape(
                R=self.min_radius,
                theta=self.end_angle * i3.DEG2RAD,
                num_points=self.num_points,
            )
            waveguide.set(
                shape=center_line,
                trace_template=self.trace_template,
                bend_radius=self.min_radius,
            )
            return waveguide

        def _generate_instances(self, insts):
            insts += i3.SRef(reference=self.waveguide)
            return insts
        
        def _generate_ports(self, ports):
            ports += self.waveguide.ports
            return ports


if __name__ == "__main__":
    cell = EulerBend()
    layout = cell.Layout(min_radius=10, end_angle=90)
    layout.visualize(annotate=True)
    layout.write_gdsii("eulerbend.gds")