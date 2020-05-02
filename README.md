[![View Skipping Stone - An interplanetary space mission design tool on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/29272-skipping-stone-an-interplanetary-space-mission-design-tool)

[![Donate to Rody](https://i.stack.imgur.com/bneea.png)](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=4M7RMVNMKAXXQ&source=url)

# Skipping Stone
A richly featured application to design trajectories for interplanetary space missions using (multiple) gravity-assist maneuvers. Its most prominent capabilities are:
- single/multi objective optimization
- can work with a multitude of (heuristic) optimization algorithms
- Low/high thrust propulsion both supported
- (rough) global and (accurate) local searches seamlessly integrated
- additional functionality can be added easily via user-definable plugins
- intuitive GUI, nice plots of best trajectory found, speed w.r.t. central body, and many others. Designed to easily explore various mission scenarios, BATCH optimization to find the best in a whole set of scenarios, etc.
- Comes ready with "nearby MPs" post processor, to also locate the known asteroids/comets/etc. that come close to the final trajectory.

I wrote this application for my Master's thesis in Space Engineering. The idea was to make everything as general and easy to understand as possible, so that other, future students could easily experiment with their problem(s), and eventually expand and improve this program. If anything, it is an excellent showcase of all that's possible with MATLAB (and all that's not :) See https://dl.dropboxusercontent.com/u/5045692/MSC_thesis_FINAL.pdf (my master thesis) for all the details and background.

To start Skipping Stone, just clone the repository (recursively, it contains submodules) or download and unpack the zip, navigate to its home directory in MATLAB and execute "main.m".

A word of advice: First run "speedup.m" (also in the home directory) first to compile a few bottleneck-algorithms; this speeds up the optimizations pretty dramatically. It is not required though; everything should work as is.

An easy example: find a so-called Mars-free-return trajectory:
- run "main.m" and leave all settings untouched
- on the "Sequence" tab, select "Venus" as swingby1, and "Mars" as target.
- Press the OPTIMIZE! button and wait a few seconds.

Feel free to add more planets and play with all the settings. Note that low thrust trajectories are MUCH harder to optimize, so you should expect to wait longer. Naturally, this also holds for increasingly complex problems; 5 swingbys with multiple revolutions and long/short way optimized will require a LOT more time than the Mars-free return mission mentioned above.

Although I took great care to catch all errors I can think of, I expect lots of bugs are still present. So if you do find any strange errors popping up, or have a much more general question, please e-mail me with a brief description, screenshots etc. so I can post a fix, or post your request here on the forum. Also, if you are also an aerospace student and want to participate, I can sure use a hand or two (or more)!

Known flaws, improvements, etc.:

- second order optimizations ("high-accuracy local search") NOT YET IMPLEMENTED: expect errors when trying this! It's something that wasn't required for my thesis, so due to time pressure at the time I didn't implement this correctly yet.

- A much better approach to building the windows and data management etc. would be to make a class out of the main window (plus all callbacks). This would make it a lot more scalable and easier to understand. But at the time I wrote this, I didn't know how to do this properly, and now that I do it's an awful lot of work: so for now this is future work.

- There is NO clear documentation (although the code is pretty well documented). This is a problem with all software written under pressure of course; it's something I'm still working on...

- As for validation: as you can read in my thesis, all missions I tried could be found by Skipping Stone, EXCEPT NASA's Galileo for some reason. Please; if you find that reason, I'd be very happy!

If you find this work useful, please consider a [donation](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=4M7RMVNMKAXXQ&source=url).
