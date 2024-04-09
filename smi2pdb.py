import os
import shutil
import time
from datetime import datetime
import static_functions
import inputs
import log_functions
from mixture import Mixture


class JobRunner:
    """
    Manages execution of computational chemistry simulations. Initializes simulation parameters,
    manages output directories, generates molecular mixtures, calculates potential energies, and logs progress.

    Attributes:
    - job_id (str): Identifier for the job, for directory and log file naming.
    - logfile (str): Path to the log file for progress and results.
    - nbr_of_mixtures_accepted (int): Count of mixtures meeting criteria.
    - trial_i (int): Current trial number.
    """

    @log_functions.track_call_depth
    def __init__(self, job_id: str, job_output_dir: str) -> None:
        """
        Initialize a JobRunner instance.

        This constructor method initializes the JobRunner instance with the provided job ID.
        It also sets initial values for instance variables, creates a log file, and prepares for the job.

        :param job_id: The identifier for the current job.
        :type job_id: str
        :return: None
        """

        self.job_id = job_id
        # Setting up job directory
        try:
            # Create the output directory for the job.
            self.job_output_dir = job_output_dir
            os.makedirs(self.job_output_dir, exist_ok=False)
        except FileExistsError:
            raise FileExistsError(
                f"Directory  already exists! Cannot submit a job there."
            )

        # Initialize instance variables with default values.
        self.logfile = self.create_log_file()
        self.nbr_of_mixtures_accepted = 0
        self.trial_i = 1

    @log_functions.track_call_depth
    def create_log_file(self) -> str:
        """
        Create a log file for the current job.

        This method creates a log file in the specified directory and returns the path
        to the created log file.

        :return: The path to the created log file.
        :rtype: str
        """
        log_file = log_functions.create_log_file(
            logdir=self.job_output_dir, constant=True
        )
        return log_file

    @log_functions.track_call_depth
    def create_directories(self) -> None:
        """
        Create directories based on the specified run mode.

        This method checks the `run_mode` specified in `inputs` and creates directories accordingly
        in the 'output' folder.

        :return: None
        """
        # Check the run mode to determine directory structure.
        os.makedirs(
            f"{self.job_output_dir}mixture_{self.nbr_of_mixtures_accepted}/pdbs/",
            exist_ok=True,
        )

    @log_functions.track_call_depth
    def run_trials(self) -> None:
        """
        Run a series of trials to create mixtures and perform simulations.

        This method iterates through multiple trials, creating mixtures and conducting simulations based
        on the specified run mode. It also logs the progress and completion time for each trial.

        :return: None
        """
        while self.trial_i <= inputs.number_of_trials:
            # Record the start time of the current trial.
            trial_start_time = time.time()

            # Create a new Mixture instance.
            mixture = Mixture()

            # Log information about the current mixture.
            log_functions.print_to_log(
                logfile=self.logfile,
                message=f"Mixture {self.nbr_of_mixtures_accepted + 1} | Initializing atomic positions (Attempt {self.trial_i})",
            )

            # Generate the mixture and proceed if accepted.
            if self.generate_mixture(mixture):
                self.nbr_of_mixtures_accepted += 1
                self.trial_i = 0
                self.create_directories()
                static_functions.write_to_pdb(
                    mol_list=mixture.molecules,
                    box_dimensions=tuple(mixture.box_dims),
                    output_file=f"{self.job_output_dir}mixture_{self.nbr_of_mixtures_accepted}/pdbs/mixture_{self.nbr_of_mixtures_accepted}.pdb",
                )

                # Exit loop if the mixtures_needed has been reached.
                if self.nbr_of_mixtures_accepted == inputs.mixtures_needed:
                    break

            # Log the completion time for the current trial.
            log_functions.print_to_log(
                logfile=self.logfile,
                message=f"Mixture {self.nbr_of_mixtures_accepted + 1} | "
                f"Mixture completed after {time.time() - trial_start_time: .2f}s.",
            )

            # Increment the trial counter.
            self.trial_i += 1

    @log_functions.track_call_depth
    def generate_mixture(self, mixture) -> bool:
        """
        Generate a mixture of molecules and calculate LJ potential energy and return True if its energy is below a threshold.

        This method generates a mixture of molecules, performs the following steps:
        1. Generates Sobol positions for the mixture.
        2. Translates molecules to Sobol positions.
        3. Calculates the Lennard-Jones potential energy for the mixture using pseudo-parameters.
        4. Returns True if the energy computed is less than 'inputs.energy_threshold'.

        :param mixture: The mixture of molecules to be generated and analyzed.
        :return: True if the mixture is accepted, False otherwise.
        :rtype: bool
        """
        log_functions.print_to_log(
            logfile=self.logfile,
            message=f"Mixture {self.nbr_of_mixtures_accepted + 1} | Generating Mixture...",
        )

        log_functions.print_to_log(
            logfile=self.logfile,
            message=f"Mixture {self.nbr_of_mixtures_accepted + 1} | Generating Sobol positions...",
        )
        # Generate Sobol positions for the mixture.
        mixture.generate_sobol_positions()

        log_functions.print_to_log(
            logfile=self.logfile,
            message=f"Mixture {self.nbr_of_mixtures_accepted + 1} | Translating molecules to Sobol positions...",
        )
        # Translate molecules to Sobol positions.
        mixture.translate_molecules_to_sobol_positions()

        log_functions.print_to_log(
            logfile=self.logfile,
            message=f"Mixture {self.nbr_of_mixtures_accepted + 1} | Calculating pseudo potential energy...",
        )
        # Calculate the LJ potential energy for the mixture.
        mixture.calculate_lj_potential_energy()

        log_functions.print_to_log(
            logfile=self.logfile,
            message=f"Mixture {self.nbr_of_mixtures_accepted + 1} | Mixture accepted? {mixture.accepted}.",
        )
        # Return whether the mixture is accepted.
        return mixture.accepted

    @log_functions.track_call_depth
    def run(self) -> None:
        """
        Run the program for creating mixtures and conducting trials.

        This method performs the following steps:
        1. Records the start time of the program.
        2. Clears the output folder associated with the job.
        3. Creates mixtures and runs trials.
        4. Logs the program's successful termination and runtime.

        :return: None
        """
        # Record the start time of the program.
        start_time = time.time()

        # Clear the output folder associated with the job.
        static_functions.clear_folder(folder_path=self.job_output_dir)

        # Print a log message indicating the start of mixture creation.
        log_functions.print_to_log(
            logfile=self.logfile, message=f"Creating mixtures..."
        )

        # Run trials to generate mixtures.
        self.run_trials()

        # Calculate and print the program's runtime.
        log_functions.print_to_log(
            logfile=self.logfile,
            message=f"Program run terminated successfully after {time.time() - start_time: .2f}s.",
        )


@log_functions.track_call_depth
def main(job_id: str, job_output_dir: str) -> None:
    """
    Main function to initialize and run a job.

    This function initializes and manages the execution of a job based on a given job ID.
    It ensures that the necessary directories are created, initializes the job runner,
    logs the initialization, and then runs the job.

    Steps:
    1. Creates the output directory for the job.
    2. Initializes a JobRunner instance for the specified job.
    3. Logs the successful initialization of the program.
    4. Runs the JobRunner to execute the job.

    To execute a simulation using only the backend, the user must call:

    Example Usage:
    ```python
    if __name__ == "__main__":
        main(job_id="1")
    ```

    The `job_id` should be an integer that uniquely identifies the job's output directory.
    The output will be stored in the '../output' directory as '../output/{job_id}/'.

    :param job_output_dir: Output directory for the job.
    :param job_id: The identifier for the current job.
    :type job_id: str
    :return: None
    :raises FileExistsError: If the output directory for the job already exists.
    """
    # Initialize a JobRunner instance for the job.
    job_runner = JobRunner(job_id, job_output_dir)

    # Log a message indicating successful program initialization.
    log_functions.print_to_log(
        logfile=job_runner.logfile, message=f"Successfully initialized program."
    )

    # Run the JobRunner to execute the job.
    job_runner.run()


if __name__ == "__main__":
    try:
        job_directory = os.getcwd()
        job_title = f"{datetime.now().strftime('%Y-%m-%d_%H:%M:%S')} | Successfully initialized program.\n"
        job_id = str(static_functions.next_job_number(directory=inputs.output_dir))
        job_output_dir = os.path.join(
            os.path.join(f"{inputs.output_dir}", ""), f"{job_id}/"
        )

        log_functions.create_function_calls_log_file(
            logfile="calls_log.txt", title=job_title
        )
        main(job_id=job_id, job_output_dir=job_output_dir)
        shutil.move("calls_log.txt", f"{job_output_dir}/calls_log.txt")

    except RuntimeError as e:
        raise RuntimeError(f"Error: {e}")
