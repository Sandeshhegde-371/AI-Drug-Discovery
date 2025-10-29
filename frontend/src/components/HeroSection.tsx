import { Beaker, Sparkles, Database, TrendingUp } from 'lucide-react';

export function HeroSection() {
  return (
    <div className="relative bg-gradient-to-br from-violet-600 via-purple-600 to-indigo-700 text-white overflow-hidden">
      {/* Background Pattern */}
      <div className="absolute inset-0 opacity-10">
        <div className="absolute inset-0" style={{
          backgroundImage: `radial-gradient(circle at 25px 25px, white 2%, transparent 0%), 
                           radial-gradient(circle at 75px 75px, white 2%, transparent 0%)`,
          backgroundSize: '100px 100px'
        }} />
      </div>
      
      {/* Floating shapes */}
      <div className="absolute top-20 left-10 w-72 h-72 bg-white/5 rounded-full blur-3xl" />
      <div className="absolute bottom-20 right-10 w-96 h-96 bg-purple-400/10 rounded-full blur-3xl" />
      
      <div className="max-w-7xl mx-auto px-6 py-20 relative z-10">
        <div className="max-w-4xl mx-auto text-center">
          {/* Badge */}
          <div className="inline-flex items-center gap-2 bg-white/10 backdrop-blur-sm border border-white/20 rounded-full px-4 py-2 mb-6">
            <Sparkles className="w-4 h-4" />
            <span className="text-sm">AI-Powered Research Platform</span>
          </div>
          
          {/* Title */}
          <h1 className="mb-6 text-5xl md:text-6xl">
            Drug Discovery Made
            <span className="block bg-gradient-to-r from-pink-200 via-purple-200 to-indigo-200 bg-clip-text text-transparent">
              Intelligent
            </span>
          </h1>
          
          {/* Subtitle */}
          <p className="text-xl opacity-90 mb-10 max-w-2xl mx-auto">
            Transform your research with AI-driven molecular generation, screening, and analysis. 
            Upload your data and unlock powerful insights in minutes.
          </p>
          
          {/* Features Grid */}
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4 max-w-3xl mx-auto">
            <div className="bg-white/10 backdrop-blur-sm border border-white/20 rounded-xl p-4">
              <Database className="w-6 h-6 mb-2 mx-auto" />
              <div className="text-sm">Flexible Data Input</div>
            </div>
            <div className="bg-white/10 backdrop-blur-sm border border-white/20 rounded-xl p-4">
              <Beaker className="w-6 h-6 mb-2 mx-auto" />
              <div className="text-sm">4 AI Pipelines</div>
            </div>
            <div className="bg-white/10 backdrop-blur-sm border border-white/20 rounded-xl p-4">
              <TrendingUp className="w-6 h-6 mb-2 mx-auto" />
              <div className="text-sm">Instant Results</div>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
