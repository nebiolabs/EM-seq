// Block comments are stripped only where `/*` opens a line. A `/*` inside a string, such as
// the glob "${params.dir}/*.fastq", would otherwise open a comment and delete every include
// up to the next `*/` anywhere in the file -- silently, leaving the hash wrong.
def stripComments(String source) {
    source.replaceAll(/(?sm)^[ \t]*\/\*.*?\*\//, '').replaceAll(/\/\/[^\n]*/, '')
}

def sha256Hex(byte[] bytes) {
    java.security.MessageDigest.getInstance('SHA-256')
        .digest(bytes)
        .collect { String.format('%02x', it) }
        .join('')
}

def resolveIncludePath(java.nio.file.Path including_file, String include_path) {
    if (include_path.contains('${')) {
        error "Include path '${include_path}' in ${including_file} is interpolated at runtime. " +
              "Include paths must be literal, because the graph is read as text, not executed."
    }

    def base = including_file.parent
    def candidates = [base.resolve(include_path), base.resolve("${include_path}.nf")]
    def resolved = candidates.find { java.nio.file.Files.isRegularFile(it) }
    if (!resolved) {
        error "Cannot resolve include '${include_path}' from ${including_file}, so the workflow content hash cannot be computed."
    }
    resolved.toRealPath()
}

def collectTransitiveIncludes(java.nio.file.Path entry_file) {
    // Matches `include { a; b as c } from 'path'`, including statements split over several lines.
    def include_pattern = ~/\binclude\s*\{[^}]*\}\s*from\s*['"]([^'"]+)['"]/
    def visited = [] as Set
    def pending = [entry_file.toRealPath()] as List

    while (pending) {
        def current_file = pending.removeLast()
        if (visited.add(current_file)) {
            def source = stripComments(current_file.getText('UTF-8'))
            (source =~ include_pattern).each { match, include_path ->
                // Anchored on current_file, not entry_file: include paths are relative to the
                // file they appear in, so anything below the first level would resolve wrongly.
                pending << resolveIncludePath(current_file, include_path)
            }
        }
    }
    visited
}

def workflowFileDigests(java.nio.file.Path entry_file, java.nio.file.Path project_dir) {
    def root = project_dir.toRealPath()
    collectTransitiveIncludes(entry_file)
        .collect { [root.relativize(it).toString(), sha256Hex(it.bytes)] }
        .sort { it[0] }
}

// Paths are hashed alongside the digests, so renaming a module changes the hash.
def workflowContentHash(java.nio.file.Path entry_file, java.nio.file.Path project_dir) {
    def digest_table = workflowFileDigests(entry_file, project_dir)
        .collect { path, digest -> "${path}\t${digest}\n" }
        .join('')
    sha256Hex(digest_table.getBytes('UTF-8'))
}

def writeWorkflowContentManifest(String workflow_name, String content_hash, java.nio.file.Path entry_file, java.nio.file.Path project_dir, String destination) {
    def digests = workflowFileDigests(entry_file, project_dir)
    def manifest = file(destination)
    manifest.parent.mkdirs()
    manifest.text = "# workflow\t${workflow_name}\n" +
                    "# content_hash\t${content_hash}\n" +
                    "# entry\t${project_dir.toRealPath().relativize(entry_file.toRealPath())}\n" +
                    "# files\t${digests.size()}\n" +
                    digests.collect { path, digest -> "${path}\t${digest}\n" }.join('')
    manifest
}
