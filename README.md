# analytical_element_model

![Preview](https://numericalenvironmental.files.wordpress.com/2016/07/particle-tracks.png)

This Python 2.7 script employs a variation on the analytic element method to simulate steady-state groundwater flow in a two-dimensional aquifer characterized by point sources/sinks, line source/sinks, faults of varying permeability, and regional flow. The steady-state model assumes a leaky aquifer is juxtaposed against an overlying aquitard of some specified thickness and vertical hydraulic conductivity. Line sources and faults are modeled using numerical integration (with a doublet element to address flow around fault segments). A particle tracking routine, constructed based upon a differential equation solver, is added as a post-processor.

The following tab-delimited input files are required (assisted with a tkinter-based GUI):

* base.txt - spatial scale and gridding (for output)
* elements.txt - definition of wells, line sources, and/or faults
* hydro.txt - hydraulic properties of aquifer and aquitard
* particlex.txt - initial positions of particles used in particle tracking routine
* track.txt - basic switches for particle tracking (on/off, forward/reverse, etc.)

More background information can be found here: https://numericalenvironmental.wordpress.com/2016/07/11/analytic-element-modeling/

Email me with questions at walt.mcnab@gmail.com. 

THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
