# References

- AI assistance (ChatGPT 5.1 Thinking, OpenAI) was used to clarify the RK4 numerical solver structure. It helped me organise the four Runge–Kutta stages and the update of the state vector (velocity, mass, angle, altitude, horizontal distance, radius) so that the integration remained consistent and numerically stable.

- AI assistance (ChatGPT 5.1 Thinking, OpenAI) supported the verification of the coupled dynamical–thermal equations. It helped me check the structure of drag, ablation, lift, gravity and curvature terms in the ODEs, and suggested ways to ensure that the implemented equations were physically consistent with the model description.

- AI assistance (ChatGPT 5.1 Thinking, OpenAI) was used to improve timestep control and stability. It provided guidance on adapting the timestep in high-gradient regions such as high-speed atmospheric entry and fragmentation onset, which reduced numerical oscillations and made the solver more robust.

- AI assistance (ChatGPT 5.1 Thinking, OpenAI) helped with coupling the pancaking (fragmentation–expansion) model into the RK4 solver. It clarified the logic for comparing dynamic pressure to material strength, refining the fragmentation trigger, and keeping the expansion of the radius and spreading rate consistent within the state updates.

- AI assistance (ChatGPT 5.1 Thinking, OpenAI) was used to refine the atmospheric-entry termination conditions. It helped formalise the criteria for classifying ground impact (altitude reaching the surface) versus complete high-altitude ablation (mass going to zero), and suggested small corrections to the final state handling just before termination.

- AI assistance (ChatGPT 5.1 Thinking, OpenAI) supported the event-type classification between airburst and cratering events. It helped organise the logic based on energy-deposition peaks to detect airbursts, and to deal with cases where no clear peak is present or where a weak event leads instead to a simple impact.

- AI assistance (ChatGPT 5.1 Thinking, OpenAI) was used to improve stability near boundary physical conditions. It provided ideas for handling extremely shallow entry angles, weak-drag regimes and situations close to the fragmentation threshold, reducing misclassification and avoiding numerical discontinuities in these edge cases.


- ChatGPT (OpenAI). Assistance was used in developing parts of the damage map generation function, including code structure and map visualisation logic. The AI-generated support included guidance on organising the function workflow, refining Folium map rendering, and
designing clear damage-zone visual elements suitable for emergency response use.

- ChatGPT (OpenAI) was used to assist in designing the probability-based circle-map visualisation, including the choice of marker styling,colour mapping, radius scaling, and folium implementation patterns.The conversation focused on improving clarity and suitability for emergency-response decision support.

- AI assistance (ChatGPT, OpenAI) was used during the development of specific components of the postcode normalisation function.
The support included clarifying the coursework-required postcode format specification (two valid 7-character formats), identifying invalid and ambiguous cases that should be rejected, helping refine the validation logic for outward and inward segments, and improving the structure and clarity of the docstring to meet documentation standards.

- AI assistance (ChatGPT, OpenAI) was used during the development of the escape-route mapping component of this project.
Specifically, the AI helped clarify how to connect OSRM routing results with Folium map layers, propose improvements to the error-handling structure for network requests, refine the logic ensuring that routes always terminate at the designated safe zone,
assist in structuring clearer and more maintainable visualisation code for risk overlays.

- AI assistance (Claude, Anthropic) was used during the optimization of the spatial search algorithm in the `find_nearest_coords` function within `locator.py`.I find out that the original KDTree-based approach could be replaced with a more efficient progressive spatial filtering method specifically tailored for uniformly distributed 1km×1km grid data.Claude assisted in implementing the bounding box filtering strategy with adaptive radius expansion, clarifying the conversion between kilometers and degrees at different latitudes, realizing efficient k-nearest neighbor selection.


- ChatGPT 5.1 Thinking (OpenAI) was used to help translate mathematical expressions into correct and readable Python code. In particular, it assisted with formulating the overpressure–distance relationships, choosing appropriate scaling laws, and embedding these equations inside helper functions such as the blast-radius calculation in a numerically robust way.

- AI assistance was used to clarify the coursework questions and marking expectations, and to explain the required tasks in more depth. It helped me break down broad instructions into smaller implementable steps, such as how `damage_zones` should use the outcome dictionary and how `impact_risk` should loop over impact scenarios and interpret the results probabilistically.

- ChatGPT 5.1 Thinking supported debugging throughout development. It helped interpret Python and Sphinx error messages, diagnose issues in function arguments and return values, and suggest safer handling of invalid or extreme inputs (e.g. zero burst energy, missing keys in the outcome, or out-of-range latitude/longitude), leading to more robust error-checking in my code.

- AI assistance was used to brainstorm edge cases and design tests for my part of the project. This included suggesting scenarios with no blast damage at the chosen pressure, runs where no postcodes or populations are returned, small synthetic outcomes to test `damage_zones`, and sampling strategies for `nsamples` in `impact_risk` so that both typical and pathological situations are exercised.

- Towards the end of implementation, ChatGPT 5.1 Thinking was used as an additional check on the overall structure and consistency of my solution. It helped me review the flow from reading the impact data, through solving atmospheric entry, computing damage zones, and aggregating postcode probabilities, and it also assisted in drafting clearer explanatory text and equations for the documentation. All final design decisions and code edits were made and verified by me.












