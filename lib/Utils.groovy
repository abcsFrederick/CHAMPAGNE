class Utils {
    // run a shell command and capture the output
    // adapted from https://groovy-lang.gitlab.io/101-scripts/basico/command_local-en.html
    public static List run(command_string) {
        def out = new StringBuilder()
        def err = new StringBuilder()
        def command = command_string.execute()
        command.consumeProcessOutput(out, err)
        command.waitFor()

        return [out, err]
    }

    // run spooker for the workflow
    public static String spooker(workflow) {
        def pipeline_name = "${workflow.manifest.name.tokenize('/')[-1]}"
        def command_string = "spooker ${workflow.launchDir} ${pipeline_name}"
        def out = new StringBuilder()
        def err = new StringBuilder()
        try {
            def command = command_string.execute()
            command.consumeProcessOutput(out, err)
            command.waitFor()
        } catch(IOError e) {
            err = e
        }
        new FileWriter("${workflow.launchDir}/log/spooker.log").with {
            write("${out}\n${err}")
            flush()
        }
        return err
    }
}
