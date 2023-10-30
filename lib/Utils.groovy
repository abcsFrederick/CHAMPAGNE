class Utils {
    // run a shell command and capture the output
    // adapted from https://groovy-lang.gitlab.io/101-scripts/basico/command_local-en.html
    public static String run(command_string) {
        def message = ""
        def out_str = new StringBuilder()
        def error_str = new StringBuilder()
        try {
            def command = command_string.execute()
            command.consumeProcessOutput(out_str, error_str)
            command.waitForOrKill(1000)

            if (!error_str.toString().equals("")){
                message = "Error executing `${command_string}`:\n${error_str}"

            } else {
                message = "Executed `${command_string}`\n${out_str}"
            }
        } catch(IOException e) {
            message = e
        }
        return message
    }

    // run spooker for the workflow
    public static String spooker(workflow) {
        def pipeline_name = "${workflow.manifest.name.tokenize('/')[-1]}"
        def message = this.run("spooker ${workflow.launchDir} ${pipeline_name} | tee ${workflow.launchDir}/log/spooker.log")
        return message
    }
}
