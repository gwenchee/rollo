import os
import subprocess
import ast
import shutil
import time
import jinja2
import time
from jinja2 import nativetypes
from rollo.openmc_evaluation import OpenMCEvaluation
from rollo.moltres_evaluation import MoltresEvaluation


class Evaluation:
    """Holds functions that generate and execute the evaluation solver's scripts.

    DEAP's (evolutionary algorithm package) fitness evaluator requires an 
    evaluation function to evaluate each individual's fitness values. The 
    Evaluation class contains a method that creates an evaluation function that 
    runs the nuclear software and returns the required fitness values, defined 
    in the input file.

    Attributes
    ----------
    supported_solvers : list of str
        list of supported evaluation software
    input_scripts : dict
        key is evaluation software name, value is that evaluation software's
        template input script name
    output_scripts : dict
        key is evaluation software name, value is that evaluation software's
        template output script name
    eval_dict : dict
        key is evaluation software name, value is a class containing the
        functions to evaluate its output files

    """

    def __init__(self):
        self.supported_solvers = ["openmc", "openmc_gc", "moltres"]
        self.input_scripts = {}
        self.output_scripts = {}
        # Developers can add new solvers to self.eval_dict below
        self.eval_dict = {"openmc": OpenMCEvaluation(), "openmc_gc": OpenMCEvaluation(), "moltres": MoltresEvaluation()}

    def add_evaluator(self, solver_name, input_script, output_script):
        """Adds information about an evaluator to the Evaluation class object
        for later use in eval_fn_generator.

        Parameters
        ----------
        solver_name : str
            name of solver
        input_script : str
            input script name
        output_script : str
            optional output script name
        """

        self.input_scripts[solver_name] = input_script
        try:
            self.output_scripts[solver_name] = output_script
        except:
            pass
        return

    def eval_fn_generator_theta(self, control_dict, output_dict, input_evaluators):
        
        def eval_fn_theta(pop):
            control_vars = self.name_ind(ind, control_dict, input_evaluators)
            order_of_solvers = self.solver_order(input_evaluators)

            all_output_vals = 1
            return all_output_vals # list of tuples
        
        return eval_fn_theta

    def eval_fn_generator(self, control_dict, output_dict, input_evaluators):
        """Returns a function that accepts a DEAP individual and returns a
        tuple of output values listed in outputs

        Parameters
        ----------
        control_dict : OrderedDict
            Ordered dict of control variables as keys and a list of their
            solver and number of variables as each value
        output_dict : OrderedDict
            Ordered dict of output variables as keys and solvers as values
        input_evaluators : dict
            evaluators sub-dictionary from input file

        Returns
        -------
        eval_function : function
            function that runs the evaluation software and returns output values
            output by the software

        """

        def eval_function(ind):
            """Accepts a DEAP individual and returns a tuple of output values
            listed in outputs

            Parameters
            ----------
            ind : deap.creator.Ind
                Created in `rollo.toolbox_generator.ToolboxGenerator`. It is
                a list with special attributes.

            Returns
            -------
            tuple
                output values from evaluators ordered by output_dict

            """

            self.rank_time = time.time()
            control_vars = self.name_ind(ind, control_dict, input_evaluators)
            output_vals = [None] * len(output_dict)
            order_of_solvers = self.solver_order(input_evaluators)
            print("ORDER", order_of_solvers)

            for solver in order_of_solvers:
                print('SOLVER',solver)
                # path name for solver's run
                path = solver + "_" + str(ind.gen) + "_" + str(ind.num)
                # render jinja-ed input script
                print(control_vars)
                if "python" in self.input_scripts[solver][0]:
                    rendered_script = self.render_jinja_template_python(
                        script=self.input_scripts[solver][1],
                        control_vars_solver=control_vars[solver],
                    )
                else:
                    rendered_script = self.render_jinja_template(
                        script=self.input_scripts[solver][1],
                        control_vars_solver=control_vars[solver],
                        ind=ind, 
                        solver=solver
                    )
                # enter directory for this particular solver's run
                os.mkdir(path)
                os.chdir(path)
                # run solver's function where run is executed
                #exec("self." + solver + "_run(rendered_script)")
                
                # run input file 
                f = open(self.input_scripts[solver][1], "w+")
                f.write(rendered_script)
                f.close()
                with open("output.txt", "wb") as output:
                    executable = self.input_scripts[solver][0].split(" ")
                    start = time.time()
                    subprocess.call(executable + [self.input_scripts[solver][1]], stdout=output)
                    end = time.time()
                    print("TIME 1",end-start)
                for i in range(len(input_evaluators[solver]["execute2"])):
                    if len(input_evaluators[solver]["execute2"][i][1]) > 1:
                        os.chdir("../")
                        shutil.copyfile(input_evaluators[solver]["execute2"][i][1], path + "/" + input_evaluators[solver]["execute2"][i][1])
                        os.chdir(path)
                        executable = input_evaluators[solver]["execute2"][i][0].split(" ") + [input_evaluators[solver]["execute2"][i][1]]
                    else:
                        executable = input_evaluators[solver]["execute2"][i][0].split(" ")

                    txt_file = "output_execute_"+str(i)+".txt"
                    with open(txt_file, "wb") as output:
                        print("execute2", i, os.getcwd())
                        start = time.time()
                        print(executable)
                        subprocess.call(executable, stdout=output)
                        end = time.time()
                        print("TIME 2",end-start)

                # go back to normal directory with all files
                os.chdir("../")
                # get output values
                if len(input_evaluators[solver]["outputs"]) > 0:
                    output_vals = self.get_output_vals(
                        output_vals, solver, output_dict, control_vars, path
                    )
                else:
                    shutil.copyfile(
                    self.output_scripts[solver][1], path + "/" + self.output_scripts[solver][1]
                )
                    # enter directory for this particular solver's run
                    os.chdir(path)
                    # run the output script
                    start = time.time()
                    oup_bytes = self.system_call(self.output_scripts[solver][0] + ' ' + self.output_scripts[solver][1])
                    end = time.time()
                    print("TIME 3",end-start)
                    # go back to normal directory with all files
                    os.chdir("../")
                if input_evaluators[solver]["keep_files"] == False:
                    shutil.rmtree(path)
            return tuple(output_vals)

        return eval_function

    def solver_order(self, input_evaluators):
        order = [None] * len(input_evaluators)
        for solver in input_evaluators:
            order[input_evaluators[solver]["order"]] = solver
        return order 

    def get_output_vals(self, output_vals, solver, output_dict, control_vars, path):
        """Returns a populated list with output values for each solver

        Parameters
        ----------
        output_vals : list
            empty list of the correct size
        solver : str
            name of solver
        output_dict : OrderedDict
            Ordered dict of output variables as keys and solvers as values
        control_vars : dict
            maps the control_dict's variable names to values from ind list to
            order the output_vals correctly
        path : str
            path name

        Returns
        -------
        output_vals : list
            output values requested by rollo input file in correct order

        """

        if self.output_scripts[solver]:
            # copy rendered output script into a new file in the particular solver's run
            shutil.copyfile(
                self.output_scripts[solver][1], path + "/" + self.output_scripts[solver][1]
            )
            # enter directory for this particular solver's run
            os.chdir(path)
            # run the output script
            start = time.time()
            oup_bytes = self.system_call(self.output_scripts[solver][0] + ' ' + self.output_scripts[solver][1])
            end = time.time()
            print("TIME 3",end-start)
            # return the output script's printed dictionary into a variable
            oup_script_results = ast.literal_eval(oup_bytes.decode("utf-8").partition('\n')[0])
            # go back to normal directory with all files
            os.chdir("../")
        for i, var in enumerate(output_dict):
            if output_dict[var] == solver:
                # if variable is a control variable
                if var in control_vars[solver]:
                    output_vals[i] = control_vars[solver][var]
                # if variable's analysis script is pre-defined
                elif var in self.eval_dict[solver].pre_defined_outputs:
                    os.chdir(path)
                    method = getattr(self.eval_dict[solver], "evaluate_" + var)
                    output_vals[i] = method()
                    os.chdir("../")
                # if variable's defined in output script
                else:
                    output_vals[i] = oup_script_results[var]
        return output_vals

    def system_call(self, command):
        """Runs evaluation software output file

        Parameters
        ----------
        command : str
            command to run in the command line

        Returns
        -------
        str
            printed output from running evaluation software output file

        """
    
        p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
        return p.stdout.read()

    def name_ind(self, ind, control_dict, input_evaluators):
        """Returns a dictionary that maps the control_dict's variable names to
        values from ind list

        Parameters
        ----------
        ind : deap.creator.Ind
            Created in `rollo.toolbox_generator.ToolboxGenerator`. It is
            a list with special attributes.
        control_dict : OrderedDict
            Ordered dict of control variables as keys and a list of their
            solver and number of variables as each value
        input_evaluators : dict
            evaluators sub-dictionary from input file

        Returns
        -------
        control_vars : dict
            maps the control_dict's variable names to values from ind list to
            order the output_vals correctly

        """

        control_vars = {}
        for solver in input_evaluators:
            control_vars[solver] = {}
        for i, var in enumerate(control_dict):
            if control_dict[var][1] == 1:
                ind_vars = ind[i]
            else:
                ind_vars = []
                for j in range(control_dict[var][1]):
                    ind_vars.append(ind[i + j])
            control_vars[control_dict[var][0]][var] = ind_vars
        return control_vars

    def render_jinja_template_python(self, script, control_vars_solver):
        """Renders a jinja2 templated python file and returns a templated python
        script. This will be used by solver's with a python interface such as
        OpenMC.

        Parameters
        ----------
        script : str
            name of evaluator template script
        control_vars_solver : str
            name of evaluation solver software

        Returns
        -------
        rendered_template : str
            rendered evaluator template script

        """

        env = nativetypes.NativeEnvironment()
        with open(script) as s:
            imported_script = s.read()
        template = nativetypes.NativeTemplate(imported_script)
        render_str = "template.render("
        for inp in control_vars_solver:
            render_str += "**{'" + inp + "':" + str(control_vars_solver[inp]) + "},"
        render_str += ")"
        rendered_template = eval(render_str)
        return rendered_template

    def render_jinja_template(self, script, control_vars_solver, ind, solver):
        """Renders a jinja2 templated input file. This will be used by solver's
        with a text based interface such as Moltres

        Parameters
        ----------
        script : str
            name of evaluator template script
        control_vars_solver : str
            name of evaluation solver software

        Returns
        -------
        rendered_template : str
            rendered evaluator template script
        """

        with open(script) as f: 
            tmp = jinja2.Template(f.read())
        render_str = "tmp.render("
        for var in control_vars_solver:
            render_str += var + "='" + control_vars_solver[var] + "',"
                # special condition for moltres 
        if solver == "moltres": 
            render_str += "group_constant_dir='../openmc_gc_" + str(ind.gen) + "_" + str(ind.num) + "'"
        render_str += ")"
        rendered_template = eval(render_str)

        return rendered_template

    def openmc_run(self, rendered_openmc_script):
        """Runs the rendered openmc script

        Parameters
        ----------
        rendered_openmc_script : str
            rendered OpenMC evaluator template script

        """

        f = open("openmc_input.py", "w+")
        f.write(rendered_openmc_script)
        f.close()
        with open("output.txt", "wb") as output:
            subprocess.call(["python", "-u", "./openmc_input.py"], stdout=output)
        return
        

    def moltres_run(self, rendered_moltres_script):
        """Runs the rendered moltres input file

        ### TO BE POPULATED
        """

        return
